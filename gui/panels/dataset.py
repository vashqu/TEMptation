"""Load data panel (blueprint Sec 8.2): one card per group, each with a
folder browser and a live pair-count preview, plus an output directory
selector. "+ Add group" is how a `treated` condition arrives later with
no code change. Calls discovery only -- never reads pixel data."""

import tkinter as tk
from tkinter import filedialog, ttk
from pathlib import Path

from temptation import discovery

from ..state import GroupEntry
from ..widgets import CARD_BG, LBL_FG, WARN_FG, make_card

DEFAULT_GROUP_NAMES = ["normal", "pathological", "treated"]


def build(parent, state):
    parent.columnconfigure(0, weight=1)

    groups_frame = tk.Frame(parent, bg=parent["bg"])
    groups_frame.grid(row=0, column=0, sticky="new")
    groups_frame.columnconfigure(0, weight=1)

    ttk.Button(
        parent, text="＋ Add group", style="Ghost.TButton",
        command=lambda: _add_group_card(groups_frame, state, _next_default_name(state)),
    ).grid(row=1, column=0, sticky="w", pady=(8, 16))

    out_outer, out_inner = make_card(parent, "Output directory")
    out_outer.grid(row=2, column=0, sticky="ew")
    out_inner.columnconfigure(1, weight=1)

    out_var = tk.StringVar(value="(not set — never defaults to the data folder)")

    def _browse_output():
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            state.output_dir = Path(path)
            out_var.set(path)
            state.mark_inputs_changed()

    tk.Label(out_inner, text="Results will be written to:", bg=CARD_BG, fg=LBL_FG).grid(
        row=0, column=0, sticky="w")
    tk.Label(out_inner, textvariable=out_var, bg=CARD_BG, fg="#999999").grid(
        row=0, column=1, sticky="w", padx=(8, 0))
    ttk.Button(out_inner, text="Browse…", style="Ghost.TButton", command=_browse_output).grid(
        row=0, column=2, sticky="e", padx=(8, 0))

    # Seed with the two conditions CLAUDE.md's data layout already has.
    _add_group_card(groups_frame, state, "normal")
    _add_group_card(groups_frame, state, "pathological")


def _next_default_name(state) -> str:
    for name in DEFAULT_GROUP_NAMES:
        if name not in state.groups:
            return name
    return f"group_{len(state.groups) + 1}"


def _add_group_card(groups_frame, state, default_name: str):
    row = len(groups_frame.grid_slaves())
    outer, inner = make_card(groups_frame, "")
    outer.grid(row=row, column=0, sticky="ew", pady=(0, 10))
    inner.columnconfigure(1, weight=1)

    entry = GroupEntry(name=default_name)
    state.groups[entry.name] = entry

    name_var = tk.StringVar(value=default_name)
    folder_var = tk.StringVar(value="(no folder selected)")
    info_var = tk.StringVar(value="")
    skip_var = tk.StringVar(value="")

    def _rename(*_a):
        new_name = name_var.get().strip() or entry.name
        if new_name == entry.name:
            return
        state.groups.pop(entry.name, None)
        entry.name = new_name
        state.groups[new_name] = entry
        state.mark_inputs_changed()

    def _browse():
        path = filedialog.askdirectory(title=f"Select folder for group '{name_var.get()}'")
        if not path:
            return
        entry.folder = Path(path)
        folder_var.set(path)
        pairs, skipped = discovery._scan_folder_recursive(Path(path))
        entry.pairs = pairs
        entry.skipped = skipped
        info_var.set(
            f"✓ {len(pairs)} pair(s) detected" if pairs else "⚠ 0 pairs detected"
        )
        skip_var.set(f"⚠ {len(skipped)} file(s) skipped  ({'; '.join(s.path.name for s in skipped[:3])}"
                      f"{', …' if len(skipped) > 3 else ''})" if skipped else "")
        state.mark_inputs_changed()

    def _remove():
        state.groups.pop(entry.name, None)
        outer.destroy()
        state.mark_inputs_changed()

    tk.Label(inner, text="Group label:", bg=CARD_BG, fg=LBL_FG).grid(row=0, column=0, sticky="w")
    name_entry = ttk.Entry(inner, textvariable=name_var, width=20)
    name_entry.grid(row=0, column=1, sticky="w", padx=(8, 0))
    name_entry.bind("<FocusOut>", _rename)

    ttk.Button(inner, text="Remove", style="Ghost.TButton", command=_remove).grid(
        row=0, column=2, sticky="e")

    tk.Label(inner, text="📁", bg=CARD_BG).grid(row=1, column=0, sticky="w", pady=(6, 0))
    tk.Label(inner, textvariable=folder_var, bg=CARD_BG, fg="#999999").grid(
        row=1, column=1, sticky="w", padx=(8, 0), pady=(6, 0))
    ttk.Button(inner, text="Browse…", style="Ghost.TButton", command=_browse).grid(
        row=1, column=2, sticky="e", pady=(6, 0))

    tk.Label(inner, textvariable=info_var, bg=CARD_BG, fg=LBL_FG,
              font=("TkDefaultFont", 9)).grid(row=2, column=0, columnspan=3, sticky="w", pady=(4, 0))
    tk.Label(inner, textvariable=skip_var, bg=CARD_BG, fg=WARN_FG,
              font=("TkDefaultFont", 9)).grid(row=3, column=0, columnspan=3, sticky="w")
