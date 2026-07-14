"""Step-rail shell for the TEMptation GUI (IMPLEMENTATION_BLUEPRINT.md
Sec 8.1). Owns the window, the rail, and the persistent status/progress
bar; delegates step content to gui/panels/*.py and all analysis work to
temptation.pipeline.analyze_dataset. This module and the panels it
imports must never compute a metric (enforced by
tests/test_gui_has_no_metrics.py, extended over this package).

Phase 7a scope (blueprint Sec 9, Phase 7 risk note: "land the rail +
dataset + calibration + run first; review/dashboard second"): steps
Metrics/QC/Review/Export are placeholders here, filled in by Phase
7b-7d without touching this file's rail/run wiring.
"""

import threading
import traceback

import tkinter as tk
from tkinter import messagebox, ttk

from temptation import pipeline

from .panels import calibration as calibration_panel
from .panels import dataset as dataset_panel
from .panels import metrics as metrics_panel
from .panels import qc as qc_panel
from .state import STEP_IDS, STEP_LABELS, AppState
from .widgets import (
    BG,
    BTN_FG,
    BTN_HOVER,
    BTN_PRIMARY,
    CARD_BG,
    ERROR_FG,
    HDR_FG,
    LBL_FG,
    STATUS_BG,
    STATUS_GLYPHS,
    apply_base_styles,
    make_card,
)

APP_TITLE = "TEMptation — Nerve Morphometry"


class TEMptationApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.configure(bg=BG)
        self.minsize(1150, 820)
        self.geometry("1300x880")
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        self.state_ = AppState()

        apply_base_styles(self)
        self._current_step = "dataset"
        self._rail_buttons: dict[str, tk.Button] = {}
        self._panel_frames: dict[str, tk.Frame] = {}

        self._build_rail()
        self._build_main_area()
        self._build_status_bar()

        self.state_.on_change(self._refresh_rail)
        self.state_.on_change(self._refresh_run_summary)

        self._select_step("dataset")
        self._refresh_rail()
        self._refresh_run_summary()

    # ------------------------------------------------------------ rail
    def _build_rail(self):
        rail = tk.Frame(self, bg="#EDEFF3", width=190)
        rail.grid(row=0, column=0, sticky="ns")
        rail.grid_propagate(False)

        tk.Label(rail, text="WORKFLOW", bg="#EDEFF3", fg="#999999",
                 font=("TkDefaultFont", 8, "bold")).pack(anchor="w", padx=16, pady=(18, 8))

        for sid in STEP_IDS:
            btn = tk.Button(
                rail, text=STEP_LABELS[sid], anchor="w",
                bg="#EDEFF3", fg=HDR_FG, activebackground="#DCE0E8",
                relief="flat", bd=0, padx=16, pady=8,
                font=("TkDefaultFont", 10),
                command=lambda s=sid: self._select_step(s),
            )
            btn.pack(fill="x")
            self._rail_buttons[sid] = btn

    def _refresh_rail(self):
        for sid, btn in self._rail_buttons.items():
            glyph, color = STATUS_GLYPHS[self.state_.step_status.get(sid, "not_started")]
            selected = sid == self._current_step
            btn.configure(
                text=f"{STEP_LABELS[sid]}   {glyph}",
                bg="#FFFFFF" if selected else "#EDEFF3",
                fg=BTN_PRIMARY if selected else HDR_FG,
            )

    def _select_step(self, step_id: str):
        self._current_step = step_id
        for sid, frame in self._panel_frames.items():
            frame.grid_remove()
        self._panel_frames[step_id].grid(row=0, column=0, sticky="nsew")
        self._refresh_rail()

    # ------------------------------------------------------------ main area
    def _build_main_area(self):
        container = tk.Frame(self, bg=BG)
        container.grid(row=0, column=1, sticky="nsew", padx=20, pady=20)
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)

        for sid in STEP_IDS:
            frame = tk.Frame(container, bg=BG)
            frame.columnconfigure(0, weight=1)
            frame.grid(row=0, column=0, sticky="nsew")
            self._panel_frames[sid] = frame

            if sid == "dataset":
                dataset_panel.build(frame, self.state_)
            elif sid == "calibration":
                calibration_panel.build(frame, self.state_)
            elif sid == "metrics":
                metrics_panel.build(frame, self.state_)
            elif sid == "qc":
                qc_panel.build(frame, self.state_)
            elif sid == "run":
                self._build_run_panel(frame)
            else:
                self._build_placeholder(frame, sid)

    def _build_placeholder(self, frame, step_id):
        tk.Label(
            frame, text=f"{STEP_LABELS[step_id]} — coming soon",
            bg=BG, fg="#999999", font=("TkDefaultFont", 12),
        ).grid(row=0, column=0, sticky="w", pady=40, padx=4)

    # ------------------------------------------------------------ run step
    def _build_run_panel(self, frame):
        outer, inner = make_card(frame, "Run analysis")
        outer.grid(row=0, column=0, sticky="ew")

        self._run_summary_var = tk.StringVar(value="")
        tk.Label(inner, textvariable=self._run_summary_var, bg=CARD_BG, fg=LBL_FG,
                 justify="left", anchor="w").pack(anchor="w", pady=(0, 12))

        self._run_btn = tk.Button(
            inner, text="▶  Run Analysis", command=self._on_run,
            bg=BTN_PRIMARY, fg=BTN_FG, activebackground=BTN_HOVER, activeforeground=BTN_FG,
            font=("TkDefaultFont", 10, "bold"), relief="flat",
            padx=14, pady=6, cursor="hand2", bd=0,
        )
        self._run_btn.pack(anchor="w")

    def _refresh_run_summary(self):
        if not hasattr(self, "_run_summary_var"):
            return
        n_groups = len(self.state_.groups)
        n_pairs = self.state_.total_pairs()
        if self.state_.pixel_space_only:
            pixel_str = "pixel-space only (no µm calibration)"
        elif self.state_.pixel_size_um:
            pixel_str = f"{self.state_.pixel_size_um:.6f} µm/px"
        else:
            pixel_str = "not set"
        out_dir = self.state_.output_dir or "(not set)"
        self._run_summary_var.set(
            f"{n_groups} group(s), {n_pairs} image pair(s) total\n"
            f"Calibration: {pixel_str}\n"
            f"Output directory: {out_dir}"
        )
        ready = (
            n_pairs > 0
            and self.state_.output_dir is not None
            and (self.state_.pixel_size_um is not None or self.state_.pixel_space_only)
            and not self.state_.run_in_progress
        )
        self._run_btn.configure(state="normal" if ready else "disabled")

    # ------------------------------------------------------------ status bar
    def _build_status_bar(self):
        bar = tk.Frame(self, bg=STATUS_BG, height=34)
        bar.grid(row=1, column=0, columnspan=2, sticky="ew")
        bar.grid_propagate(False)
        bar.columnconfigure(0, weight=1)

        self._progress = ttk.Progressbar(bar, mode="determinate", maximum=100)
        self._progress.grid(row=0, column=0, sticky="ew", padx=(12, 8), pady=8)

        self._status_var = tk.StringVar(value="Ready")
        tk.Label(bar, textvariable=self._status_var, bg=STATUS_BG, fg=LBL_FG,
                 font=("TkDefaultFont", 9)).grid(row=0, column=1, padx=8)

        self._cancel_btn = tk.Button(
            bar, text="Cancel", command=self._on_cancel,
            bg=STATUS_BG, fg=ERROR_FG, relief="flat", bd=0, state="disabled",
        )
        self._cancel_btn.grid(row=0, column=2, padx=(4, 12))

    # ------------------------------------------------------------ run / cancel
    def _on_run(self):
        pairs = self.state_.all_pairs()
        if not pairs:
            messagebox.showwarning("No data", "Add at least one group with detected pairs first.")
            return
        if self.state_.output_dir is None:
            messagebox.showwarning("No output directory",
                                    "Set an output directory in the Load data step first.")
            return
        if self.state_.pixel_size_um is None and not self.state_.pixel_space_only:
            messagebox.showwarning("No calibration",
                                    "Set a pixel size in the Calibrate step, "
                                    "or tick 'Report pixel-space metrics only'.")
            return

        self.state_.run_in_progress = True
        self.state_.cancel_event = threading.Event()
        self._run_btn.configure(state="disabled")
        self._cancel_btn.configure(state="normal")
        self._progress.configure(value=0, maximum=len(pairs))
        self._status_var.set(f"Processing 0/{len(pairs)}…")

        # Pixel-space-only mode isn't yet a first-class backend schema
        # variant (blueprint Sec 8.3 / CLAUDE.md Sec 16.3 call for
        # _px-suffixed, _um2-free output -- future work); until then this
        # runs with an identity scale and the checkbox is a documented
        # placeholder, not silently misleading math.
        pixel_um = self.state_.pixel_size_um or 1.0

        threading.Thread(
            target=self._run_worker,
            args=(pairs, pixel_um, self.state_.cancel_event),
            daemon=True,
        ).start()

    def _on_cancel(self):
        self.state_.cancel_event.set()
        self._status_var.set("Cancelling…")

    def _run_worker(self, pairs, pixel_um, cancel_event):
        def on_done(i, total, pair, *_rest):
            self.after(0, lambda: self._on_progress(i + 1, total, pair["id"]))

        def on_error(i, total, pair, exc, tb_str):
            self.after(0, lambda: self._on_progress(i + 1, total, pair["id"], error=str(exc)))

        try:
            result = pipeline.analyze_dataset(
                pairs,
                pixel_length_um=pixel_um,
                seg_cfg=self.state_.seg_cfg,
                qc_thresholds=self.state_.qc_thresholds,
                exclude_qc_failed=self.state_.exclude_qc_failed,
                collect_mito=self.state_.collect_mito,
                collect_qc_report=True,
                on_image_done=on_done,
                on_image_error=on_error,
                cancel_event=cancel_event,
            )
            self.after(0, lambda: self._run_done(result, cancel_event))
        except Exception:
            tb = traceback.format_exc()
            self.after(0, lambda: self._run_failed(tb))

    def _on_progress(self, done, total, image_id, error=None):
        self._progress.configure(value=done, maximum=total)
        if error:
            self._status_var.set(f"Processing {done}/{total}: {image_id} — ERROR: {error}")
        else:
            self._status_var.set(f"Processing {done}/{total}: {image_id}")

    def _run_done(self, result, cancel_event):
        self.state_.run_in_progress = False
        self.state_.result = result
        self.state_.last_raw_axons = result.df_axons
        self._cancel_btn.configure(state="disabled")
        n_err = len(result.errors)
        if cancel_event.is_set():
            self._status_var.set(f"Cancelled — {len(result.df_axons)} axon(s) from completed images kept")
            self.state_.step_status["run"] = "warning"
        elif n_err:
            self._status_var.set(
                f"Done with {n_err} error(s) — {len(result.df_axons)} axons, {len(result.df_image)} images")
            self.state_.step_status["run"] = "warning"
        else:
            self._status_var.set(f"Done — {len(result.df_axons)} axons, {len(result.df_image)} images")
            self.state_.step_status["run"] = "complete"
        self.state_.step_status["review"] = "ready"
        self.state_.step_status["export"] = "ready"
        self.state_.notify()
        self._refresh_run_summary()

    def _run_failed(self, tb: str):
        self.state_.run_in_progress = False
        self._cancel_btn.configure(state="disabled")
        self._status_var.set("Run failed — see dialog")
        self._refresh_run_summary()
        messagebox.showerror("Analysis failed", f"An error occurred:\n\n{tb}")


def main():
    app = TEMptationApp()
    app.mainloop()


if __name__ == "__main__":
    main()
