"""Metric selection panel (blueprint Sec 8.4): seven informational cards
plus an advanced accordion of SegmentationConfig fields, reflected from
the dataclass so a new field there gets a GUI row with no GUI edit.

The backend always computes every metric a loaded mask's compartments
support -- there is no per-metric compute toggle to wire up. What this
panel actually controls is (a) whether mitochondria_metrics.csv gets
collected during Run (state.collect_mito, the one real "advanced
mitochondrial metrics" switch CLAUDE.md Sec 16.4 asks for) and (b)
mask-availability gating, so a card for metrics a group's masks can't
support is visibly explained rather than just silently full of NaN.

Mask scanning is an explicit "Scan mask availability" button rather
than something re-run on every state change -- re-reading every mask
in a 100-image batch on each calibration keystroke elsewhere in the
app would make unrelated steps janky."""

import dataclasses
import tkinter as tk
from tkinter import ttk

from temptation import dataio
from temptation.config import SegmentationConfig

from ..widgets import CARD_BG, HDR_FG, LBL_FG, SUCCESS_FG, WARN_FG, Tooltip, make_card

# (title, required mask kinds, one-line rationale, always computed)
METRIC_CARDS = [
    ("Core axon & myelin", (), "Axon and fiber area, perimeter, circularity.", True),
    ("g-ratio distribution", ("myelin",), "CV / median / IQR / min / max of g-ratio per image.", False),
    ("Mitochondrial burden", ("mito",), "Occupancy, density, fragmentation, normalized load.", False),
    ("Mitochondrial shape", ("mito",), "Aspect ratio, solidity, eccentricity (skimage regionprops).", False),
    ("Spatial distribution", ("mito",), "Peripheralization index and nearest-neighbor clustering.", False),
    ("Axon-myelin relationship", ("myelin", "mito"), "normalized_mito_load, mito_per_myelin.", False),
    ("Image-level summaries", (), "Writes image_summary.csv with per-image aggregates.", True),
]


def build(parent, state):
    parent.columnconfigure(0, weight=1)

    tk.Label(
        parent,
        text=("Every metric below is always computed for any image whose mask has "
              "the required compartment -- this panel doesn't turn computation on "
              "or off, it tells you which metrics will be meaningful for your "
              "loaded data, and which will be NaN."),
        bg=parent["bg"], fg="#999999", font=("TkDefaultFont", 9),
        wraplength=640, justify="left",
    ).grid(row=0, column=0, sticky="w", pady=(0, 10))

    collect_mito_var = tk.BooleanVar(value=state.collect_mito)

    def _on_collect_mito_toggle():
        state.collect_mito = collect_mito_var.get()
        state.mark_inputs_changed()

    ttk.Checkbutton(
        parent, text="Collect mitochondria_metrics.csv (one row per mitochondrion) during Run",
        style="Card.TCheckbutton", variable=collect_mito_var, command=_on_collect_mito_toggle,
    ).grid(row=1, column=0, sticky="w", pady=(0, 12))

    scan_row = tk.Frame(parent, bg=parent["bg"])
    scan_row.grid(row=2, column=0, sticky="w", pady=(0, 12))
    availability_var = tk.StringVar(value="Mask availability not checked yet.")
    ttk.Button(scan_row, text="Scan mask availability", style="Ghost.TButton",
               command=lambda: _scan(state, availability_var, card_widgets)).pack(side="left")
    tk.Label(scan_row, textvariable=availability_var, bg=parent["bg"], fg=LBL_FG,
             font=("TkDefaultFont", 9)).pack(side="left", padx=(10, 0))

    cards_frame = tk.Frame(parent, bg=parent["bg"])
    cards_frame.grid(row=3, column=0, sticky="ew")
    cards_frame.columnconfigure((0, 1), weight=1)

    card_widgets = {}
    for i, (title, required, rationale, always_on) in enumerate(METRIC_CARDS):
        outer, inner = make_card(cards_frame, title)
        outer.grid(row=i // 2, column=i % 2, sticky="ew",
                   padx=(0, 10) if i % 2 == 0 else (0, 0), pady=(0, 10))
        badge = "requires: " + ", ".join(required) if required else "no extra mask required"
        tk.Label(inner, text=badge, bg=CARD_BG, fg="#999999", font=("TkDefaultFont", 8)).pack(anchor="w")
        status_var = tk.StringVar(value="always on" if always_on else "not yet scanned")
        status_lbl = tk.Label(inner, textvariable=status_var, bg=CARD_BG, fg=HDR_FG, font=("TkDefaultFont", 9))
        status_lbl.pack(anchor="w", pady=(4, 0))
        Tooltip(outer, rationale)
        card_widgets[title] = (status_var, status_lbl, required, always_on)

    _build_advanced_accordion(parent, state)


def _scan(state, availability_var, card_widgets):
    pairs = state.all_pairs()
    if not pairs:
        availability_var.set("No masks loaded yet -- add a group in Load data first.")
        return

    n = len(pairs)
    n_myelin = 0
    n_mito = 0
    for p in pairs:
        try:
            mask = dataio.read_mask(p["mask_path"])
        except Exception:
            continue
        hist = dataio.mask_value_histogram(mask)
        if hist.get(state.seg_cfg.myelin_val, 0) > 0:
            n_myelin += 1
        if hist.get(state.seg_cfg.mito_val, 0) > 0:
            n_mito += 1

    availability_var.set(f"Scanned {n} mask(s): myelin in {n_myelin}, mitochondria in {n_mito}.")

    presence = {"myelin": n_myelin, "mito": n_mito}
    for title, (status_var, status_lbl, required, always_on) in card_widgets.items():
        if always_on:
            continue
        missing = [kind for kind in required if presence[kind] == 0]
        partial = [kind for kind in required if 0 < presence[kind] < n]
        if missing:
            status_var.set("⚠ no " + " or ".join(missing) + " pixels in any loaded mask")
            status_lbl.configure(fg=WARN_FG)
        elif partial:
            status_var.set(
                "; ".join(f"{kind} present in {presence[kind]}/{n} images" for kind in partial)
            )
            status_lbl.configure(fg=WARN_FG)
        else:
            status_var.set("available for all loaded images")
            status_lbl.configure(fg=SUCCESS_FG)


def _build_advanced_accordion(parent, state):
    outer, inner = make_card(parent, "Advanced segmentation parameters")
    outer.grid(row=4, column=0, sticky="ew", pady=(14, 0))
    inner.columnconfigure(1, weight=1)

    fields = dataclasses.fields(SegmentationConfig)
    entries = {}
    for i, f in enumerate(fields):
        current = getattr(state.seg_cfg, f.name)
        var = tk.StringVar(value=str(current))
        tk.Label(inner, text=f.name, bg=CARD_BG, fg=LBL_FG, font=("TkDefaultFont", 9)).grid(
            row=i, column=0, sticky="w", pady=1)
        entry = ttk.Entry(inner, textvariable=var, width=12)
        entry.grid(row=i, column=1, sticky="w", padx=(8, 0), pady=1)
        entries[f.name] = (var, type(current))

    error_var = tk.StringVar(value="")

    def _apply():
        kwargs = {}
        bad = []
        for name, (var, orig_type) in entries.items():
            text = var.get().strip()
            try:
                if orig_type is bool:
                    kwargs[name] = text.lower() in ("1", "true", "yes")
                else:
                    kwargs[name] = orig_type(text)
            except ValueError:
                kwargs[name] = getattr(state.seg_cfg, name)
                bad.append(name)
        state.seg_cfg = SegmentationConfig(**kwargs)
        error_var.set(f"Kept previous value for: {', '.join(bad)}" if bad else "")
        state.mark_inputs_changed()

    ttk.Button(inner, text="Apply", style="Ghost.TButton", command=_apply).grid(
        row=len(fields), column=0, columnspan=2, sticky="w", pady=(8, 0))
    tk.Label(inner, textvariable=error_var, bg=CARD_BG, fg=WARN_FG, font=("TkDefaultFont", 8)).grid(
        row=len(fields) + 1, column=0, columnspan=2, sticky="w")
