"""QC sub-tab of the "Setup" step (blueprint Sec 8.5): threshold fields
reflected from config.QCThresholds (a new threshold field there gets a
GUI row with no GUI edit), plus a live preview of flag counts
recomputed via temptation.qc.compute_qc_flags -- against the last run's
already-computed raw per-axon metrics (state.result.df_axons, which
survives later parameter edits -- see gui/state.py's
mark_inputs_changed), never by re-segmenting and never by
reimplementing a flag rule here (that recomputation is the entire
reason "qc" is an approved temptation submodule for this panel).

Threshold edits are debounced: committing state.qc_thresholds (and so
calling mark_inputs_changed()/notify()) on every keystroke used to
cascade into every listener across the app -- rebuilding the dashboard
figure, the export checklist, etc. -- dozens of times while typing a
single number."""

import dataclasses
import tkinter as tk
from tkinter import ttk

from temptation import qc as qc_mod
from temptation.config import QCThresholds

from ..widgets import CARD_BG, LBL_FG, debounce, make_card


def build(parent, state):
    parent.columnconfigure(0, weight=1)

    outer, inner = make_card(parent, "QC thresholds")
    outer.grid(row=0, column=0, sticky="ew")
    inner.columnconfigure(1, weight=1)

    fields = dataclasses.fields(QCThresholds)
    entries = {}
    for i, f in enumerate(fields):
        current = getattr(state.qc_thresholds, f.name)
        var = tk.StringVar(value="" if current is None else str(current))
        tk.Label(inner, text=f.name, bg=CARD_BG, fg=LBL_FG, font=("TkDefaultFont", 9)).grid(
            row=i, column=0, sticky="w", pady=2)
        entry = ttk.Entry(inner, textvariable=var, width=12)
        entry.grid(row=i, column=1, sticky="w", padx=(8, 0), pady=2)
        numeric_type = int if isinstance(current, int) and not isinstance(current, bool) else float
        entries[f.name] = (var, numeric_type)

    exclude_var = tk.BooleanVar(value=state.exclude_qc_failed)

    preview_var = tk.StringVar(value="")

    def _apply(*_a):
        kwargs = {}
        for name, (var, numeric_type) in entries.items():
            text = var.get().strip()
            if text == "":
                kwargs[name] = None
                continue
            try:
                kwargs[name] = numeric_type(text)
            except ValueError:
                kwargs[name] = getattr(state.qc_thresholds, name)
        state.qc_thresholds = QCThresholds(**kwargs)
        state.exclude_qc_failed = exclude_var.get()
        state.mark_inputs_changed()
        _refresh_preview()

    debounced_apply = debounce(inner, 400, _apply)
    for name, (var, _t) in entries.items():
        var.trace_add("write", lambda *_a: debounced_apply())

    cb = ttk.Checkbutton(
        inner, text="Exclude QC-failed objects from image summaries",
        style="Card.TCheckbutton", variable=exclude_var, command=_apply,
    )
    cb.grid(row=len(fields), column=0, columnspan=2, sticky="w", pady=(10, 0))
    tk.Label(
        inner, text="Raw values are always exported. Excluded objects are marked, never deleted.",
        bg=CARD_BG, fg="#999999", font=("TkDefaultFont", 8),
    ).grid(row=len(fields) + 1, column=0, columnspan=2, sticky="w")

    preview_outer, preview_inner = make_card(parent, "Live preview")
    preview_outer.grid(row=1, column=0, sticky="ew", pady=(14, 0))
    tk.Label(preview_inner, textvariable=preview_var, bg=CARD_BG, fg=LBL_FG,
             justify="left", anchor="w", font=("TkDefaultFont", 9)).pack(anchor="w")

    last_seen = {"key": None}

    def _refresh_preview():
        result = state.result
        df = result.df_axons if result is not None else None
        if df is None or df.empty:
            preview_var.set("Run analysis once (① Setup) to see a live QC preview here.")
            last_seen["key"] = None
            return

        # QCThresholds is a frozen dataclass (value-equal, not just
        # identity-equal) -- gate on (result, thresholds) so an
        # unrelated notify() elsewhere in Setup (adding a Data group,
        # toggling collect_mito) doesn't re-run compute_qc_flags over
        # every axon for no reason. _apply() calls this directly right
        # after committing a real threshold change, so genuine edits
        # still refresh immediately.
        key = (id(result), state.qc_thresholds)
        if key == last_seen["key"]:
            return
        last_seen["key"] = key

        flagged = qc_mod.compute_qc_flags(df, state.qc_thresholds)
        qc_cols = [c for c in flagged.columns if c.startswith("qc_")]
        n_images = df["image_id"].nunique() if "image_id" in df.columns else 1
        lines = [f"{len(df)} axon(s) across {n_images} image(s)"]
        any_flagged = False
        for col in qc_cols:
            n = int(flagged[col].sum())
            if n:
                any_flagged = True
                pct = 100 * n / len(df)
                lines.append(f"  {col}: {n} ({pct:.0f}%)")
        if not any_flagged:
            lines.append("  no flags triggered at current thresholds")
        preview_var.set("\n".join(lines))

    state.on_change(_refresh_preview)
    _refresh_preview()
