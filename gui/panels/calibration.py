"""Calibration sub-tab of the "Setup" step (blueprint Sec 8.3): a single
pixel-size field with a derived readout, and the pixel-space-only
escape hatch (CLAUDE.md Sec 16.3 -- the user must not be able to run
calibrated metrics without either a pixel size or an explicit opt-out).
The `1/value` readout below is a display convenience the user typed
the input for, not a derived biological metric, so it stays out of the
temptation package.

The derived readout/state.pixel_size_um update on every keystroke
(cheap, local), but the state.mark_inputs_changed() call -- which
fans out to every panel's state.on_change listener -- is debounced, so
typing a pixel size doesn't repeatedly demote the rail and rebuild the
dashboard/export checklist mid-keystroke."""

import tkinter as tk
from tkinter import ttk

from ..widgets import CARD_BG, ERROR_FG, LBL_FG, Tooltip, debounce, make_card


def build(parent, state):
    parent.columnconfigure(0, weight=1)

    outer, inner = make_card(parent, "Calibration")
    outer.grid(row=0, column=0, sticky="ew")
    inner.columnconfigure(1, weight=1)

    pixel_var = tk.StringVar(value="" if state.pixel_size_um is None else str(state.pixel_size_um))
    derived_var = tk.StringVar(value="")
    error_var = tk.StringVar(value="")
    pixel_only_var = tk.BooleanVar(value=state.pixel_space_only)

    lbl = tk.Label(inner, text="Pixel size (µm / pixel):", bg=CARD_BG, fg=LBL_FG)
    lbl.grid(row=0, column=0, sticky="w")
    entry = ttk.Entry(inner, textvariable=pixel_var, width=14)
    entry.grid(row=0, column=1, sticky="w", padx=(8, 0))
    tip = ("All area, distance, and density metrics scale with this value.\n"
           "g-ratio, circularity, and solidity do not.\n"
           "Derive from a scale bar: pixel_size_um = bar_length_um / bar_length_px.")
    Tooltip(lbl, tip)
    Tooltip(entry, tip)

    tk.Label(inner, textvariable=derived_var, bg=CARD_BG, fg="#999999",
             font=("TkDefaultFont", 9)).grid(row=1, column=0, columnspan=2, sticky="w", pady=(4, 0))
    tk.Label(inner, textvariable=error_var, bg=CARD_BG, fg=ERROR_FG,
             font=("TkDefaultFont", 9)).grid(row=2, column=0, columnspan=2, sticky="w")

    debounced_notify = debounce(entry, 400, state.mark_inputs_changed)

    def _on_pixel_change(*_a):
        text = pixel_var.get().strip()
        if not text:
            state.pixel_size_um = None
            derived_var.set("")
            error_var.set("")
        else:
            try:
                value = float(text)
                if not value > 0:
                    raise ValueError
            except ValueError:
                state.pixel_size_um = None
                derived_var.set("")
                error_var.set("Pixel size must be a positive number.")
            else:
                state.pixel_size_um = value
                error_var.set("")
                derived_var.set(f"1 µm ≈ {1 / value:.1f} px")
        debounced_notify()

    entry.bind("<KeyRelease>", _on_pixel_change)
    _on_pixel_change()

    def _on_pixel_only_toggle():
        state.pixel_space_only = pixel_only_var.get()
        entry.configure(state="disabled" if state.pixel_space_only else "normal")
        state.mark_inputs_changed()

    cb = ttk.Checkbutton(
        inner, text="Report pixel-space metrics only (no calibration available)",
        style="Card.TCheckbutton", variable=pixel_only_var, command=_on_pixel_only_toggle,
    )
    cb.grid(row=3, column=0, columnspan=2, sticky="w", pady=(10, 0))
    if state.pixel_space_only:
        entry.configure(state="disabled")

    tk.Label(
        inner,
        text="All area, distance, and density metrics scale with this value; "
             "g-ratio, circularity, and solidity do not.",
        bg=CARD_BG, fg="#999999", font=("TkDefaultFont", 8),
        wraplength=520, justify="left",
    ).grid(row=4, column=0, columnspan=2, sticky="w", pady=(10, 0))
