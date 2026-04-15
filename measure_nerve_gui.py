#!/usr/bin/env python3

"""
measure_nerve_gui.py
--------------------
Tkinter GUI for nerve fiber segmentation & measurement.

Standalone companion to the CLI script `measure_nerve.py`.
Both share the same core function: `measure_image()`.

Requirements (same as CLI + tkinter which ships with Python):
    numpy, pandas, scipy, scikit-image, tifffile, matplotlib

Usage:
    python measure_nerve_gui.py
"""

import sys
import threading
import traceback
from pathlib import Path

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Import core engine from the CLI script (must be on PYTHONPATH or same dir)
# ---------------------------------------------------------------------------
_this_dir = Path(__file__).resolve().parent
if str(_this_dir) not in sys.path:
    sys.path.insert(0, str(_this_dir))

from measure_nerve import measure_image, _find_pairs_in_folder  # noqa: E402


# ---------------------------------------------------------------------------
# Design tokens
# ---------------------------------------------------------------------------
BG          = "#F5F7FA"
CARD_BG     = "white"
HDR_FG      = "#333333"
LBL_FG      = "#666666"
BORDER_CLR  = "#CCCCCC"
BTN_PRIMARY = "#007BFF"
BTN_HOVER   = "#0056b3"
BTN_FG      = "white"
STATUS_BG   = "#EFEFEF"
EMPTY_FG    = "#888888"
ROW_ODD     = "white"
ROW_EVEN    = "#F9F9F9"

APP_TITLE   = "Nerve Fiber Analyzer"
CARD_PAD    = 14


# ---------------------------------------------------------------------------
# Tooltip helper
# ---------------------------------------------------------------------------

class Tooltip:
    """Show a small floating tooltip when the mouse hovers over a widget."""

    def __init__(self, widget: tk.Widget, text: str):
        self._widget = widget
        self._text   = text
        self._tip: tk.Toplevel | None = None
        widget.bind("<Enter>", self._show, add="+")
        widget.bind("<Leave>", self._hide, add="+")

    def _show(self, _event=None):
        if self._tip or not self._text:
            return
        x = self._widget.winfo_rootx() + 20
        y = self._widget.winfo_rooty() + self._widget.winfo_height() + 4
        self._tip = tw = tk.Toplevel(self._widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        tk.Label(
            tw,
            text=self._text,
            justify="left",
            background="#ffffc0",
            foreground="#222222",
            relief="solid",
            borderwidth=1,
            font=("TkDefaultFont", 9),
            wraplength=340,
            padx=6,
            pady=5,
        ).pack()

    def _hide(self, _event=None):
        if self._tip:
            self._tip.destroy()
            self._tip = None


# ---------------------------------------------------------------------------
# Card helper
# ---------------------------------------------------------------------------

def _make_card(parent, title: str, **grid_kw) -> tuple[tk.Frame, tk.Frame]:
    """White card with a bold title.  Returns (outer_frame, inner_frame)."""
    outer = tk.Frame(parent, bg=CARD_BG,
                     highlightthickness=1, highlightbackground=BORDER_CLR)
    outer.grid(**grid_kw)
    tk.Label(outer, text=title, bg=CARD_BG, fg=HDR_FG,
             font=("TkDefaultFont", 11, "bold")).pack(
        anchor="w", padx=CARD_PAD, pady=(CARD_PAD, 4))
    inner = tk.Frame(outer, bg=CARD_BG)
    inner.pack(fill="both", expand=True, padx=CARD_PAD, pady=(0, CARD_PAD))
    return outer, inner


# ---------------------------------------------------------------------------
# Main application
# ---------------------------------------------------------------------------

class NerveApp(tk.Tk):
    """Main Tkinter application window."""

    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.configure(bg=BG)
        self.minsize(1100, 850)
        self.geometry("1200x900")
        self.resizable(True, True)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)   # status bar

        # State
        self._tem_path: str | None      = None
        self._mask_path: str | None     = None
        self._folder_path: str | None   = None
        self._batch_pairs: list         = []
        self._df_axons: pd.DataFrame | None  = None
        self._df_image: pd.DataFrame | None  = None
        self._labels_ws: np.ndarray | None   = None
        self._resolved_mode: str | None      = None
        self._tem_array: np.ndarray | None   = None
        self._mask_array: np.ndarray | None  = None
        # Must be created after super().__init__() and before _build_ui()
        self._batch_mode_var = tk.StringVar(value="single")

        self._apply_styles()
        self._build_ui()
        self._build_status_bar()

    # ------------------------------------------------------------------ styles
    def _apply_styles(self):
        s = ttk.Style(self)
        try:
            s.theme_use("clam")
        except Exception:
            pass

        s.configure("TFrame",            background=BG)
        s.configure("Card.TFrame",       background=CARD_BG)

        s.configure("TLabel",            background=BG,      foreground=LBL_FG)
        s.configure("Card.TLabel",       background=CARD_BG, foreground=LBL_FG)

        s.configure("TEntry",
                    fieldbackground="white",
                    bordercolor=BORDER_CLR,
                    lightcolor=BORDER_CLR,
                    darkcolor=BORDER_CLR,
                    insertcolor=HDR_FG)

        s.configure("Ghost.TButton",
                    background=CARD_BG,
                    foreground="#444444",
                    borderwidth=1,
                    relief="solid",
                    padding=(8, 4))
        s.map("Ghost.TButton",
              background=[("active", "#F0F0F0"), ("disabled", CARD_BG)],
              foreground=[("disabled", "#AAAAAA")])

        s.configure("TCombobox",
                    fieldbackground="white",
                    bordercolor=BORDER_CLR)

        s.configure("TNotebook",         background=BG, tabmargins=[2, 5, 2, 0])
        s.configure("TNotebook.Tab",
                    background="#DDE1E7",
                    foreground=HDR_FG,
                    padding=[12, 5])
        s.map("TNotebook.Tab",
              background=[("selected", CARD_BG)],
              foreground=[("selected", BTN_PRIMARY)])

        s.configure("Treeview",
                    background=ROW_ODD,
                    fieldbackground=ROW_ODD,
                    rowheight=22,
                    foreground=HDR_FG)
        s.configure("Treeview.Heading",
                    background="#EAECEF",
                    foreground=HDR_FG,
                    relief="flat",
                    font=("TkDefaultFont", 9, "bold"))
        s.map("Treeview",
              background=[("selected", BTN_PRIMARY)],
              foreground=[("selected", "white")])

        s.configure("Card.TCheckbutton", background=CARD_BG, foreground=LBL_FG)

    # ------------------------------------------------------------------ status bar
    def _build_status_bar(self):
        bar = tk.Frame(self, bg=STATUS_BG, height=26)
        bar.grid(row=1, column=0, sticky="ew")
        bar.grid_propagate(False)
        self._status_var = tk.StringVar(value="Ready")
        tk.Label(bar, textvariable=self._status_var,
                 bg=STATUS_BG, fg="#666666",
                 font=("TkDefaultFont", 9)).pack(side="right", padx=10, pady=4)

    # ------------------------------------------------------------------ main UI
    def _build_ui(self):
        main = tk.Frame(self, bg=BG)
        main.grid(row=0, column=0, sticky="nsew", padx=15, pady=15)
        main.columnconfigure(0, weight=1)
        main.rowconfigure(3, weight=1)

        self._build_file_card(main)
        self._build_param_card(main)
        self._build_action_row(main)
        self._build_notebook(main)

    # ------------------------------------------------------------------ file card
    def _build_file_card(self, parent):
        _, inner = _make_card(parent, "Image Files",
                              row=0, column=0, sticky="ew", pady=(0, 12))
        inner.columnconfigure(1, weight=1)

        # ---- Mode selector (row=0) ----
        mode_frame = tk.Frame(inner, bg=CARD_BG)
        mode_frame.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        tk.Label(mode_frame, text="Input mode:", bg=CARD_BG, fg=LBL_FG).pack(
            side="left", padx=(0, 8))
        tk.Radiobutton(
            mode_frame, text="Single image", bg=CARD_BG, fg=LBL_FG,
            activebackground=CARD_BG, selectcolor=CARD_BG,
            variable=self._batch_mode_var, value="single",
            command=self._toggle_mode,
        ).pack(side="left")
        tk.Radiobutton(
            mode_frame, text="Batch folder", bg=CARD_BG, fg=LBL_FG,
            activebackground=CARD_BG, selectcolor=CARD_BG,
            variable=self._batch_mode_var, value="batch",
            command=self._toggle_mode,
        ).pack(side="left", padx=(10, 0))

        # ---- Single-image frame (row=1) ----
        self._single_frame = tk.Frame(inner, bg=CARD_BG)
        self._single_frame.columnconfigure(1, weight=1)

        # TEM row
        tk.Label(self._single_frame, text="TEM image:", bg=CARD_BG, fg=LBL_FG).grid(
            row=0, column=0, sticky="w", pady=(0, 4))
        self._tem_var = tk.StringVar(value="(none — optional, needed for overlay plot)")
        tk.Label(self._single_frame, textvariable=self._tem_var,
                 bg=CARD_BG, fg="#999999").grid(
            row=0, column=1, sticky="w", padx=(8, 0), pady=(0, 4))
        btn_t = tk.Frame(self._single_frame, bg=CARD_BG)
        btn_t.grid(row=0, column=2, sticky="e", padx=(8, 0), pady=(0, 4))
        ttk.Button(btn_t, text="Browse…", style="Ghost.TButton",
                   command=self._browse_tem).pack(side="left")
        ttk.Button(btn_t, text="Clear",   style="Ghost.TButton",
                   command=self._clear_tem).pack(side="left", padx=(4, 0))

        # Mask row
        tk.Label(self._single_frame, text="Mask image:", bg=CARD_BG, fg=LBL_FG).grid(
            row=1, column=0, sticky="w")
        self._mask_var = tk.StringVar(value="(none)")
        tk.Label(self._single_frame, textvariable=self._mask_var,
                 bg=CARD_BG, fg="#999999").grid(
            row=1, column=1, sticky="w", padx=(8, 0))
        btn_m = tk.Frame(self._single_frame, bg=CARD_BG)
        btn_m.grid(row=1, column=2, sticky="e", padx=(8, 0))
        ttk.Button(btn_m, text="Browse…", style="Ghost.TButton",
                   command=self._browse_mask).pack(side="left")

        # Mask info row (unique values, shown after a mask is loaded)
        self._mask_info_var = tk.StringVar(value="")
        tk.Label(self._single_frame, textvariable=self._mask_info_var,
                 bg=CARD_BG, fg="#888888",
                 font=("TkDefaultFont", 8),
                 anchor="w", justify="left").grid(
            row=2, column=0, columnspan=3, sticky="ew", pady=(2, 0))

        # ---- Batch folder frame (row=1, same slot — one hidden at a time) ----
        self._batch_frame = tk.Frame(inner, bg=CARD_BG)
        self._batch_frame.columnconfigure(1, weight=1)

        # Folder row
        tk.Label(self._batch_frame, text="Folder:", bg=CARD_BG, fg=LBL_FG).grid(
            row=0, column=0, sticky="w", pady=(0, 4))
        self._folder_var = tk.StringVar(value="(none — select folder containing TEM + mask pairs)")
        tk.Label(self._batch_frame, textvariable=self._folder_var,
                 bg=CARD_BG, fg="#999999").grid(
            row=0, column=1, sticky="w", padx=(8, 0), pady=(0, 4))
        btn_f = tk.Frame(self._batch_frame, bg=CARD_BG)
        btn_f.grid(row=0, column=2, sticky="e", padx=(8, 0), pady=(0, 4))
        ttk.Button(btn_f, text="Browse folder…", style="Ghost.TButton",
                   command=self._browse_folder).pack(side="left")

        # Folder info row
        self._folder_info_var = tk.StringVar(value="")
        tk.Label(self._batch_frame, textvariable=self._folder_info_var,
                 bg=CARD_BG, fg="#888888",
                 font=("TkDefaultFont", 8),
                 anchor="w", justify="left").grid(
            row=1, column=0, columnspan=3, sticky="ew", pady=(2, 0))

        # Grid both content frames into row=1; hide batch frame by default
        self._single_frame.grid(row=1, column=0, columnspan=3, sticky="ew")
        self._batch_frame.grid(row=1, column=0, columnspan=3, sticky="ew")
        self._batch_frame.grid_remove()

    # ------------------------------------------------------------------ mode toggle
    def _toggle_mode(self):
        if self._batch_mode_var.get() == "batch":
            self._single_frame.grid_remove()
            self._batch_frame.grid()
            self._plot_btn.configure(state="disabled")
        else:
            self._batch_frame.grid_remove()
            self._single_frame.grid()

    # ------------------------------------------------------------------ param card
    def _build_param_card(self, parent):
        _, inner = _make_card(parent, "Parameters",
                              row=1, column=0, sticky="ew", pady=(0, 12))

        def _param(pf, row, text, var, tip="", col=0, w=10):
            lbl = tk.Label(pf, text=text, bg=CARD_BG, fg=LBL_FG)
            lbl.grid(row=row, column=col, sticky="w", padx=(0, 2))
            ent = ttk.Entry(pf, textvariable=var, width=w)
            ent.grid(row=row, column=col + 1, sticky="w", padx=(2, 10), pady=3)
            if tip:
                Tooltip(lbl, tip)
                Tooltip(ent, tip)
            return ent

        def _combo(pf, row, text, var, vals, tip="", col=0, w=10):
            lbl = tk.Label(pf, text=text, bg=CARD_BG, fg=LBL_FG)
            lbl.grid(row=row, column=col, sticky="w", padx=(0, 2))
            cb = ttk.Combobox(pf, textvariable=var, values=vals,
                              state="readonly", width=w)
            cb.grid(row=row, column=col + 1, sticky="w", padx=(2, 10), pady=3)
            if tip:
                Tooltip(lbl, tip)
                Tooltip(cb, tip)
            return cb

        # ---- Row 0: pixel size | mode | watershed | ws weight | ws beta ----
        self._pixel_um_var = tk.StringVar(value="0.00524")
        _param(inner, 0, "Pixel size (µm/px):", self._pixel_um_var,
               tip=(
                   "Pixel size in micrometers per pixel.\n"
                   "Default: 0.00524\n"
                   "Derived from scale bar: divide scale bar length (µm) by its width in pixels.\n"
                   "Example: bar = 191 px → 1 µm  ⟹  1 / 191 ≈ 0.00524 µm/px"
               ), col=0)

        self._mode_var = tk.StringVar(value="auto")
        _combo(inner, 0, "Mode:", self._mode_var,
               ["auto", "normal", "pathological"],
               tip=(
                   "Segmentation mode.\n"
                   "  auto        — auto-detects normal vs pathological based on myelin pixel count\n"
                   "  normal      — expects myelin sheath; computes g-ratio, myelin thickness, MVF\n"
                   "  pathological — no myelin expected; axon-only measurements\n"
                   "Default: auto"
               ), col=2, w=12)

        self._ws_mode_var = tk.StringVar(value="weighted")
        _combo(inner, 0, "Watershed:", self._ws_mode_var,
               ["weighted", "simple"],
               tip=(
                   "Watershed algorithm variant for assigning myelin area to axons.\n"
                   "  weighted — larger axons expand further into shared territory\n"
                   "  simple   — standard distance-transform watershed (equal expansion)\n"
                   "Default: weighted"
               ), col=4, w=10)

        self._ws_weight_var = tk.StringVar(value="radius")
        _combo(inner, 0, "WS weight:", self._ws_weight_var,
               ["radius", "area"],
               tip=(
                   "Weighting metric for the weighted watershed.\n"
                   "  radius — weights by equivalent radius (equivalent_diameter / 2)\n"
                   "  area   — weights by √(axon area)\n"
                   "Default: radius"
               ), col=6, w=8)

        self._ws_beta_var = tk.StringVar(value="1.0")
        _param(inner, 0, "Watershed β (beta):", self._ws_beta_var,
               tip=(
                   "β controls the bias favoring larger axons during weighted watershed.\n"
                   "Elevation is computed as: elevation = max(distance − β·weight, 0).\n"
                   "Higher β makes large axons claim more myelin territory; "
                   "too high β can over-bias boundaries."
               ), col=8, w=5)

        # ---- Row 1: mask label values + smoothing ----
        self._myelin_val_var   = tk.StringVar(value="64")
        self._axoplasm_val_var = tk.StringVar(value="192")
        self._mito_val_var     = tk.StringVar(value="128")
        self._smoothing_var    = tk.StringVar(value="1")

        _param(inner, 1, "Myelin val:", self._myelin_val_var,
               tip=(
                   "Pixel intensity value that represents myelin in the mask image.\n"
                   "Default: 64\n"
                   "Must match the value used when the mask was created."
               ), col=0)
        _param(inner, 1, "Axoplasm val:", self._axoplasm_val_var,
               tip=(
                   "Pixel intensity value that represents axoplasm in the mask image.\n"
                   "Default: 192\n"
                   "Must match the value used when the mask was created."
               ), col=2)
        _param(inner, 1, "Mito val:", self._mito_val_var,
               tip=(
                   "Pixel intensity value that represents mitochondria in the mask image.\n"
                   "Default: 128\n"
                   "Must match the value used when the mask was created."
               ), col=4)
        _param(inner, 1, "Smoothing (px):", self._smoothing_var,
               tip=(
                   "Morphological smoothing radius in pixels applied to binary masks\n"
                   "(binary opening on axoplasm, binary closing on myelin).\n"
                   "Default: 1\n"
                   "Range: 0 (no smoothing) – 5\n"
                   "Set to 0 to disable smoothing entirely."
               ), col=6)

        # ---- Row 2: area thresholds + compactness ----
        self._min_axon_var      = tk.StringVar(value="200")
        self._min_myelin_var    = tk.StringVar(value="300")
        self._myelin_thresh_var = tk.StringVar(value="200")
        self._compactness_var   = tk.StringVar(value="0.001")

        _param(inner, 2, "Min axon area (px):", self._min_axon_var,
               tip=(
                   "Minimum axon area in pixels. Objects smaller than this are\n"
                   "discarded as noise before watershed.\n"
                   "Default: 200\n"
                   "Increase if spurious tiny objects appear; decrease if small axons are missed."
               ), col=0)
        _param(inner, 2, "Min myelin area (px):", self._min_myelin_var,
               tip=(
                   "Minimum myelin area in pixels. Below this threshold,\n"
                   "myelin-dependent metrics (g-ratio, thickness) are set to NaN.\n"
                   "The axon row is still kept in the output.\n"
                   "Default: 300"
               ), col=2)
        _param(inner, 2, "Myelin threshold (px):", self._myelin_thresh_var,
               tip=(
                   "Pixel count threshold used in 'auto' mode to decide the segmentation mode.\n"
                   "If total myelin pixel count ≥ threshold → normal mode; otherwise → pathological.\n"
                   "Default: 200"
               ), col=4)
        _param(inner, 2, "Compactness:", self._compactness_var,
               tip=(
                   "Watershed compactness parameter. Higher values produce rounder,\n"
                   "more regular watershed regions (penalises irregular boundaries).\n"
                   "Default: 0.001\n"
                   "Range: 0 (natural boundaries) – 1.0 (very compact/circular)"
               ), col=6)

        # ---- Row 3: detached myelin checkbox ----
        self._assign_myelin_var = tk.BooleanVar(value=False)
        cb = ttk.Checkbutton(
            inner,
            text="Assign Detached Myelin to Nearest Axon",
            variable=self._assign_myelin_var,
            style="Card.TCheckbutton",
        )
        cb.grid(row=3, column=0, columnspan=8, sticky="w", pady=(8, 0))
        Tooltip(cb,
                "When enabled, myelin pixels outside the watershed fiber regions\n"
                "are assigned to the nearest axon by Euclidean distance.\n"
                "Allows computing myelin metrics for axons with detached sheaths.\n"
                "Default: off (use connected/watershed myelin only)")

    # ------------------------------------------------------------------ action row
    def _build_action_row(self, parent):
        row = tk.Frame(parent, bg=BG)
        row.grid(row=2, column=0, sticky="ew", pady=(0, 12))

        # Primary button — tk.Button for colour support
        self._run_btn = tk.Button(
            row,
            text="▶  Run Analysis",
            command=self._on_run,
            bg=BTN_PRIMARY,
            fg=BTN_FG,
            activebackground=BTN_HOVER,
            activeforeground=BTN_FG,
            font=("TkDefaultFont", 10, "bold"),
            relief="flat",
            padx=14,
            pady=6,
            cursor="hand2",
            bd=0,
        )
        self._run_btn.pack(side="left")
        self._run_btn.bind(
            "<Enter>",
            lambda _e: self._run_btn.configure(bg=BTN_HOVER)
            if str(self._run_btn["state"]) != "disabled" else None,
        )
        self._run_btn.bind(
            "<Leave>",
            lambda _e: self._run_btn.configure(bg=BTN_PRIMARY)
            if str(self._run_btn["state"]) != "disabled" else None,
        )

        # Ghost secondary buttons
        self._plot_btn = ttk.Button(row, text="Show Overlay Plot",
                                    style="Ghost.TButton",
                                    command=self._on_plot, state="disabled")
        self._plot_btn.pack(side="left", padx=(8, 0))

        self._export_btn = ttk.Button(row, text="Export CSV…",
                                      style="Ghost.TButton",
                                      command=self._on_export, state="disabled")
        self._export_btn.pack(side="left", padx=(8, 0))

    # ------------------------------------------------------------------ notebook
    def _build_notebook(self, parent):
        nb = ttk.Notebook(parent)
        nb.grid(row=3, column=0, sticky="nsew")

        axon_tab = tk.Frame(nb, bg=CARD_BG)
        nb.add(axon_tab, text="Axon Metrics")
        axon_tab.columnconfigure(0, weight=1)
        axon_tab.rowconfigure(0, weight=1)
        self._axon_tree = self._make_table(axon_tab)
        self._axon_tree.grid(row=0, column=0, sticky="nsew")

        summary_tab = tk.Frame(nb, bg=CARD_BG)
        nb.add(summary_tab, text="Image Summary")
        summary_tab.columnconfigure(0, weight=1)
        summary_tab.rowconfigure(0, weight=1)
        self._summary_tree = self._make_table(summary_tab)
        self._summary_tree.grid(row=0, column=0, sticky="nsew")

    # ------------------------------------------------------------------ table
    @staticmethod
    def _make_table(parent) -> tk.Frame:
        outer = tk.Frame(parent, bg=CARD_BG)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(0, weight=1)

        # Tree sub-frame (shown when data is present)
        tree_frame = tk.Frame(outer, bg=CARD_BG)
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)

        tree = ttk.Treeview(tree_frame, show="headings", selectmode="browse")
        vsb  = ttk.Scrollbar(tree_frame, orient="vertical",   command=tree.yview)
        hsb  = ttk.Scrollbar(tree_frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.tag_configure("oddrow",  background=ROW_ODD)
        tree.tag_configure("evenrow", background=ROW_EVEN)

        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        # Empty-state frame (shown when no data)
        empty_frame = tk.Frame(outer, bg=CARD_BG)
        tk.Label(empty_frame, text="⬡", fg="#CCCCCC", bg=CARD_BG,
                 font=("TkDefaultFont", 36)).pack(pady=(60, 10))
        tk.Label(
            empty_frame,
            text=(
                "Ready to Analyze.\n"
                "Upload a TEM image and mask, then click\n"
                "Run Analysis to view metrics here."
            ),
            fg=EMPTY_FG, bg=CARD_BG,
            font=("TkDefaultFont", 10),
            justify="center",
        ).pack()

        # Start with empty state visible
        empty_frame.grid(row=0, column=0, sticky="nsew")

        # Stash references as attributes on the container
        outer._tree        = tree        # type: ignore[attr-defined]
        outer._tree_frame  = tree_frame  # type: ignore[attr-defined]
        outer._empty_frame = empty_frame # type: ignore[attr-defined]

        return outer

    def _populate_table(self, outer: tk.Frame, df: pd.DataFrame):
        tree        = outer._tree        # type: ignore[attr-defined]
        tree_frame  = outer._tree_frame  # type: ignore[attr-defined]
        empty_frame = outer._empty_frame # type: ignore[attr-defined]

        tree.delete(*tree.get_children())

        if df is None or df.empty:
            tree_frame.grid_remove()
            empty_frame.grid(row=0, column=0, sticky="nsew")
            return

        empty_frame.grid_remove()
        tree_frame.grid(row=0, column=0, sticky="nsew")

        cols = list(df.columns)
        tree["columns"] = cols

        MIN_COL_W = 70
        MAX_COL_W = 150
        CHAR_PX   = 7

        for c in cols:
            w = max(MIN_COL_W, min(MAX_COL_W, len(c) * CHAR_PX + 16))
            tree.heading(c, text=c, anchor="w")
            tree.column(c, width=w, minwidth=MIN_COL_W, stretch=False)

        for i, (_, row) in enumerate(df.iterrows()):
            vals = []
            for v in row:
                if isinstance(v, float):
                    vals.append(f"{v:.4g}" if not np.isnan(v) else "—")
                else:
                    vals.append(str(v))
            tag = "evenrow" if i % 2 == 0 else "oddrow"
            tree.insert("", "end", values=vals, tags=(tag,))

    # ------------------------------------------------------------------ browsing
    _FILETYPES = [
        ("TIFF images", "*.tif *.tiff"),
        ("PNG images",  "*.png"),
        ("All files",   "*.*"),
    ]

    def _browse_tem(self):
        path = filedialog.askopenfilename(title="Select TEM image",
                                          filetypes=self._FILETYPES)
        if path:
            self._tem_path = path
            self._tem_var.set(Path(path).name)

    def _clear_tem(self):
        self._tem_path = None
        self._tem_var.set("(none — optional, needed for overlay plot)")

    def _browse_mask(self):
        path = filedialog.askopenfilename(title="Select Mask image",
                                          filetypes=self._FILETYPES)
        if path:
            self._mask_path = path
            self._mask_var.set(Path(path).name)
            self._update_mask_info(path)

    def _update_mask_info(self, path: str):
        """Read the mask file and display unique pixel values with counts."""
        try:
            import tifffile as tiff
            arr = tiff.imread(path)
        except Exception:
            try:
                from PIL import Image
                arr = np.array(Image.open(path))
            except Exception:
                self._mask_info_var.set("(could not read mask for value info)")
                return
        if arr.ndim > 2:
            arr = arr[..., 0]
        unique, counts = np.unique(arr, return_counts=True)
        parts = [f"{int(v)}: {int(c):,} px" for v, c in zip(unique, counts)]
        self._mask_info_var.set("Detected values:  " + "   ".join(parts))

    def _browse_folder(self):
        path = filedialog.askdirectory(title="Select folder with TEM + mask image pairs")
        if path:
            self._folder_path = path
            self._folder_var.set(Path(path).name)
            self._update_folder_info(path)

    def _update_folder_info(self, path: str):
        """Scan folder for TEM+mask pairs and show a summary."""
        try:
            pairs = _find_pairs_in_folder(Path(path))
            self._batch_pairs = pairs
            if pairs:
                ids_preview = ", ".join(p["id"] for p in pairs[:6])
                suffix = "…" if len(pairs) > 6 else ""
                self._folder_info_var.set(
                    f"Found {len(pairs)} TEM+mask pair(s): {ids_preview}{suffix}"
                )
            else:
                self._folder_info_var.set(
                    "No TEM+mask pairs found. "
                    "Files must contain 'axon'/'tem' or 'mask' in their names."
                )
        except Exception as e:
            self._folder_info_var.set(f"Error scanning folder: {e}")

    # ------------------------------------------------------------------ params
    def _read_params(self) -> dict:
        """Read all GUI fields into a dict.  Raises ValueError on bad input."""
        return {
            "pixel_length_um":        float(self._pixel_um_var.get()),
            "myelin_val":             int(self._myelin_val_var.get()),
            "axoplasm_val":           int(self._axoplasm_val_var.get()),
            "mito_val":               int(self._mito_val_var.get()),
            "mode":                   self._mode_var.get(),
            "myelin_threshold_px":    int(self._myelin_thresh_var.get()),
            "smoothing_radius_px":    int(self._smoothing_var.get()),
            "min_axon_area_px":       int(self._min_axon_var.get()),
            "min_myelin_area_px":     int(self._min_myelin_var.get()),
            "watershed_mode":         self._ws_mode_var.get(),
            "watershed_weight":       self._ws_weight_var.get(),
            "watershed_compactness":  float(self._compactness_var.get()),
            "watershed_beta":         float(self._ws_beta_var.get()),
            "assign_detached_myelin": "nearest" if self._assign_myelin_var.get() else "none",
        }

    # ------------------------------------------------------------------ run
    def _on_run(self):
        if self._batch_mode_var.get() == "batch":
            if not self._folder_path:
                messagebox.showwarning("No folder", "Please select a folder first.")
                return
            if not self._batch_pairs:
                messagebox.showwarning(
                    "No pairs",
                    "No TEM+mask pairs found in the selected folder.\n"
                    "Files must contain 'axon'/'tem' or 'mask' in their names.",
                )
                return
        else:
            if not self._mask_path:
                messagebox.showwarning("No mask", "Please select a mask image first.")
                return

        try:
            params = self._read_params()
        except ValueError as e:
            messagebox.showerror("Bad parameter", f"Invalid parameter value:\n{e}")
            return

        self._run_btn.configure(state="disabled", bg="#6CAED8")
        self._status_var.set("Running…")
        self.update_idletasks()

        if self._batch_mode_var.get() == "batch":
            threading.Thread(
                target=self._run_batch_analysis, args=(params,), daemon=True
            ).start()
        else:
            threading.Thread(
                target=self._run_analysis, args=(params,), daemon=True
            ).start()

    def _run_analysis(self, params: dict):
        try:
            import tifffile as tiff

            mask = tiff.imread(self._mask_path)
            self._mask_array = mask

            if self._tem_path:
                tem = tiff.imread(self._tem_path)
                self._tem_array = tem
            else:
                tem = np.zeros(mask.shape[:2], dtype=np.uint8)
                self._tem_array = None

            df_axons, df_image, labels_ws, resolved_mode = measure_image(
                tem=tem, mask=mask, **params
            )

            self._df_axons      = df_axons
            self._df_image      = df_image
            self._labels_ws     = labels_ws
            self._resolved_mode = resolved_mode

            self.after(0, self._analysis_done)

        except Exception:
            tb = traceback.format_exc()
            self.after(0, lambda: self._analysis_error(tb))

    def _analysis_done(self):
        n        = len(self._df_axons) if self._df_axons is not None else 0
        mode_str = self._resolved_mode or "?"
        self._status_var.set(f"Done — {n} axons detected  [{mode_str}]")
        self._run_btn.configure(state="normal", bg=BTN_PRIMARY)
        self._export_btn.configure(state="normal" if n > 0 else "disabled")

        if self._tem_array is not None and n > 0:
            self._plot_btn.configure(state="normal")
        else:
            self._plot_btn.configure(state="disabled")
            if self._tem_array is None and n > 0:
                messagebox.showinfo(
                    "Metrics only",
                    "No TEM image was loaded, so only metrics are available.\n"
                    "Load a TEM image and re-run to also get the overlay plot."
                )

        self._populate_table(self._axon_tree,    self._df_axons)
        self._populate_table(self._summary_tree, self._df_image)

    def _analysis_error(self, tb: str):
        self._status_var.set("Error — see dialog")
        self._run_btn.configure(state="normal", bg=BTN_PRIMARY)
        messagebox.showerror("Analysis failed", f"An error occurred:\n\n{tb}")

    # ------------------------------------------------------------------ batch
    def _run_batch_analysis(self, params: dict):
        try:
            import tifffile as tiff

            pairs = _find_pairs_in_folder(Path(self._folder_path))
            n = len(pairs)
            all_axons  = []
            all_images = []

            for i, pair in enumerate(pairs):
                img_id = pair["id"]
                self.after(
                    0,
                    lambda msg=f"Processing {i + 1}/{n}: {img_id}…": (
                        self._status_var.set(msg)
                    ),
                )

                tem_arr  = tiff.imread(pair["tem_path"])
                mask_arr = tiff.imread(pair["mask_path"])

                df_axons, df_image, _, resolved_mode = measure_image(
                    tem=tem_arr, mask=mask_arr, **params
                )

                for df in (df_axons, df_image):
                    if not df.empty:
                        if "image_id" not in df.columns:
                            df.insert(0, "image_id", img_id)
                        if "mode" not in df.columns:
                            df.insert(1, "mode", resolved_mode)

                if not df_axons.empty:
                    all_axons.append(df_axons)
                if not df_image.empty:
                    all_images.append(df_image)

            self._df_axons = (
                pd.concat(all_axons, ignore_index=True) if all_axons else pd.DataFrame()
            )
            self._df_image = (
                pd.concat(all_images, ignore_index=True) if all_images else pd.DataFrame()
            )
            self._labels_ws     = None   # no single label image in batch mode
            self._resolved_mode = "batch"
            self._tem_array     = None
            self._mask_array    = None

            total = len(self._df_axons)
            self.after(0, lambda: self._batch_done(n, total))

        except Exception:
            tb = traceback.format_exc()
            self.after(0, lambda: self._analysis_error(tb))

    def _batch_done(self, n_images: int, total_axons: int):
        self._status_var.set(
            f"Done — {n_images} image(s) processed, {total_axons} axons total"
        )
        self._run_btn.configure(state="normal", bg=BTN_PRIMARY)
        self._export_btn.configure(state="normal" if total_axons > 0 else "disabled")
        self._plot_btn.configure(state="disabled")  # no single overlay in batch mode
        self._populate_table(self._axon_tree,    self._df_axons)
        self._populate_table(self._summary_tree, self._df_image)

    # ------------------------------------------------------------------ plot
    def _on_plot(self):
        if self._tem_array is None:
            messagebox.showinfo("No TEM image",
                                "Overlay plot requires a TEM image.\n"
                                "Load a TEM and re-run to enable plotting.")
            return

        if self._labels_ws is None or self._df_axons is None or self._df_axons.empty:
            messagebox.showinfo("No data", "Run analysis first.")
            return

        try:
            params = self._read_params()
        except ValueError:
            params = {}

        try:
            import matplotlib
            matplotlib.use("TkAgg")
            import matplotlib.pyplot as plt
            from scipy import ndimage as ndi
            from skimage.morphology import binary_dilation, disk
            from skimage import measure as sk_measure
            from matplotlib.patches import Patch
            from matplotlib.lines import Line2D

            tem       = self._tem_array
            mask      = self._mask_array
            labels_ws = self._labels_ws
            df_axons  = self._df_axons

            myelin_val   = params.get("myelin_val",   64)
            axoplasm_val = params.get("axoplasm_val", 192)
            mito_val     = params.get("mito_val",     128)

            MYELIN_RGBA   = np.array([0.27, 0.51, 0.71, 0.45], dtype=np.float32)
            AXOPLASM_RGBA = np.array([1.00, 0.60, 0.10, 0.50], dtype=np.float32)
            MITO_RGBA     = np.array([0.85, 0.15, 0.15, 0.80], dtype=np.float32)

            fig, axes = plt.subplots(1, 2, figsize=(16, 7))
            fig.suptitle(f"Nerve Fiber Analysis  [{self._resolved_mode}]", fontsize=13)

            # Build 2-D arrays for pixel-hover info (format_coord)
            _mask_2d = mask if mask.ndim == 2 else mask[..., 0]
            _tem_2d  = tem  if tem.ndim  == 2 else tem[..., 0]
            _h, _w   = _mask_2d.shape

            def _make_fmt(include_tem: bool):
                def _fmt(x, y):
                    ix, iy = int(round(x)), int(round(y))
                    if 0 <= iy < _h and 0 <= ix < _w:
                        mv  = int(_mask_2d[iy, ix])
                        msg = f"x={ix}  y={iy}  mask={mv}"
                        if include_tem:
                            tv  = int(_tem_2d[iy, ix])
                            msg += f"  tem={tv}"
                        return msg
                    return f"x={x:.0f}  y={y:.0f}"
                return _fmt

            axes[0].imshow(tem, cmap="gray", interpolation="nearest")
            axes[0].format_coord = _make_fmt(include_tem=True)
            axes[0].set_title("TEM image")
            axes[0].axis("off")

            axes[1].imshow(tem, cmap="gray", interpolation="nearest")
            axes[1].format_coord = _make_fmt(include_tem=True)

            overlay = np.zeros((*labels_ws.shape, 4), dtype=np.float32)
            # Tissue overlays based solely on raw mask values, NOT restricted by
            # watershed labels_ws.  labels_ws is used only for contours and axon IDs.
            axon_only     = (mask == axoplasm_val)
            mito_raw      = (mask == mito_val)
            axon_dil      = binary_dilation(axon_only, disk(1))
            mito_in_axon  = mito_raw & axon_dil
            axoplasm_full = axon_only | mito_in_axon

            myelin_mask   = (mask == myelin_val)
            axoplasm_mask = axoplasm_full.copy()
            mito_mask     = ndi.binary_fill_holes(mito_raw)
            axoplasm_mask = axoplasm_mask & ~mito_mask

            overlay[myelin_mask]   = MYELIN_RGBA
            overlay[axoplasm_mask] = AXOPLASM_RGBA
            overlay[mito_mask]     = MITO_RGBA
            axes[1].imshow(overlay, interpolation="nearest")

            unique_labels = np.unique(labels_ws)
            unique_labels = unique_labels[unique_labels > 0]
            for lbl in unique_labels:
                fiber_bin = (labels_ws == lbl).astype(np.uint8)
                contours  = sk_measure.find_contours(fiber_bin, level=0.5)
                for contour in contours:
                    axes[1].plot(contour[:, 1], contour[:, 0],
                                 color="white", linewidth=0.8, alpha=0.9)

            if not df_axons.empty:
                for _, row in df_axons.iterrows():
                    axes[1].text(
                        row["centroid_x_px"], row["centroid_y_px"],
                        str(int(row["axon_id"])),
                        color="white", fontsize=6, ha="center", va="center",
                        fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.45, lw=0),
                    )

            legend_elements = [
                Patch(facecolor=MYELIN_RGBA[:3],   alpha=0.8, label="Myelin"),
                Patch(facecolor=AXOPLASM_RGBA[:3], alpha=0.8, label="Axoplasm"),
                Patch(facecolor=MITO_RGBA[:3],     alpha=0.8, label="Mitochondria"),
                Line2D([0], [0], color="white", linewidth=1.2, label="Fiber boundary"),
            ]
            axes[1].legend(handles=legend_elements, loc="lower right",
                           fontsize=7, framealpha=0.6)
            axes[1].set_title("Overlay + axon IDs")
            axes[1].axis("off")
            plt.tight_layout()
            plt.show(block=True)

        except Exception:
            tb = traceback.format_exc()
            messagebox.showerror("Plot error", f"Could not create plot:\n\n{tb}")

    # ------------------------------------------------------------------ export
    def _on_export(self):
        if self._df_axons is None or self._df_axons.empty:
            messagebox.showinfo("No data", "Run analysis first.")
            return

        path = filedialog.asksaveasfilename(
            title="Save axon metrics CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile="axons.csv",
        )
        if not path:
            return

        try:
            self._df_axons.to_csv(path, index=False)
            summary_path = Path(path).with_name(Path(path).stem + "_image_summary.csv")
            if self._df_image is not None and not self._df_image.empty:
                self._df_image.to_csv(summary_path, index=False)
                messagebox.showinfo(
                    "Exported",
                    f"Axon metrics → {path}\nImage summary → {summary_path}"
                )
            else:
                messagebox.showinfo("Exported", f"Axon metrics → {path}")
        except Exception as e:
            messagebox.showerror("Export failed", str(e))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    app = NerveApp()
    app.mainloop()


if __name__ == "__main__":
    main()
