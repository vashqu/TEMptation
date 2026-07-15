"""Visual review panel (blueprint Sec 8.6) -- the direct answer to
segmentation artifacts otherwise looking like biology (the class of bug
that motivated Phase 2b's mitochondrial-hole fix in the first place).

Re-analyzes exactly one selected image on demand via
pipeline.analyze_image_legacy: analyze_dataset doesn't retain every
image's labels_ws for a whole batch (a deliberate memory tradeoff from
the Phase 7 prep work), so this panel re-derives it for just the image
being looked at, which is cheap (well under a second per image).

This file draws nothing itself -- it decides which already-computed
layer temptation.plotting.overlay_figure should draw and what color an
already-computed axon_id's boundary should be, and does simple
nearest-centroid click hit-testing (generic UI geometry, not a
biological formula). That is why plotting.overlay_figure gained
show_*/boundary_colors/df_mito parameters in this same phase, instead
of this panel calling skimage.measure.find_contours itself.

Layout note: overlay_figure()'s default figsize (16, 7) is sized for
its own popup window, not this 3-column embed (list | canvas |
inspector). Passed straight through, the canvas's natural width
crowded the list/inspector columns down to near-zero once an image
loaded -- the actual cause of "the image list disappears after I pick
one" (the matplotlib toolbar's Home/Back/Forward buttons some users
then reach for are pan/zoom-history controls, unrelated to switching
images -- clicking them doing nothing is expected once the list is
gone, not a separate bug). Fixed here with a smaller figsize plus
minsize on the outer columns; Previous/Next buttons give a second,
layout-independent way to change images."""

import tkinter as tk
from tkinter import ttk

import numpy as np

from temptation import dataio, pipeline

from ..widgets import CARD_BG, ERROR_FG, LBL_FG, make_card

QC_BADGE_GLYPH = {"normal": "●", "flagged": "⚠", "excluded": "✕", "unknown": "○"}
QC_BOUNDARY_COLOR = {"normal": "#AAAAAA", "flagged": "#F6AD55", "excluded": ERROR_FG}


def build(parent, state):
    # minsize on the outer columns is the actual fix for the "image list
    # disappears" bug: without it, the center canvas's natural width
    # (matplotlib's default figsize) squeezed columns 0/2 toward zero.
    parent.columnconfigure(0, minsize=190)
    parent.columnconfigure(1, weight=1)
    parent.columnconfigure(2, minsize=230)
    parent.rowconfigure(0, weight=1)

    current = {"pair": None, "df_axons": None, "df_mito": None,
               "tem": None, "mask": None, "labels_ws": None}
    canvas_holder = {"canvas": None, "toolbar": None, "fig": None}

    # ---------------------------------------------------------- left: image list
    list_outer, list_inner = make_card(parent, "Images")
    list_outer.grid(row=0, column=0, sticky="ns", padx=(0, 12))
    listbox = tk.Listbox(list_inner, width=30, height=30, exportselection=False)
    listbox.pack(fill="both", expand=True)

    # ---------------------------------------------------------- center: overlay
    center = tk.Frame(parent, bg=parent["bg"])
    center.grid(row=0, column=1, sticky="nsew")
    center.rowconfigure(1, weight=1)
    center.columnconfigure(0, weight=1)

    toggles_row = tk.Frame(center, bg=parent["bg"])
    toggles_row.grid(row=0, column=0, sticky="w", pady=(0, 6))

    nav_row = tk.Frame(center, bg=parent["bg"])
    nav_row.grid(row=0, column=0, sticky="e", pady=(0, 6))

    layer_vars = {
        "show_myelin": tk.BooleanVar(value=True),
        "show_axoplasm": tk.BooleanVar(value=True),
        "show_mito": tk.BooleanVar(value=True),
        "show_boundaries": tk.BooleanVar(value=True),
        "show_axon_ids": tk.BooleanVar(value=True),
    }
    qc_color_mode = tk.BooleanVar(value=True)
    show_mito_sublayer = tk.BooleanVar(value=True)

    canvas_frame = tk.Frame(center, bg=CARD_BG)
    canvas_frame.grid(row=1, column=0, sticky="nsew")
    canvas_frame.rowconfigure(0, weight=1)
    canvas_frame.columnconfigure(0, weight=1)

    # ---------------------------------------------------------- right: inspector
    inspector_outer, inspector_inner = make_card(parent, "Axon inspector")
    inspector_outer.grid(row=0, column=2, sticky="ns", padx=(12, 0))
    inspector_var = tk.StringVar(value="Select an image, then click an axon in the overlay to inspect it.")
    tk.Label(inspector_inner, textvariable=inspector_var, bg=CARD_BG, fg=LBL_FG,
             justify="left", anchor="w", wraplength=230, font=("TkDefaultFont", 9)).pack(anchor="w")

    # ---------------------------------------------------------- helpers

    def _image_qc_status(image_id) -> str:
        """Checks both axon-level qc_* flags (df_axons) and image-level
        ones (df_image, e.g. qc_mask_noncanonical/F6) -- an image can be
        worth a second look even when every individual axon in it looks
        fine, so a badge sourced from df_axons alone would miss it."""
        result = state.result
        if result is None or result.df_axons.empty or "image_id" not in result.df_axons.columns:
            return "unknown"
        sub = result.df_axons[result.df_axons["image_id"].astype(str) == str(image_id)]
        if sub.empty:
            return "unknown"
        if "excluded_from_analysis" in sub.columns and sub["excluded_from_analysis"].any():
            return "excluded"

        axon_qc_cols = [c for c in sub.columns if c.startswith("qc_")]
        flagged = bool(axon_qc_cols and sub[axon_qc_cols].any(axis=1).any())

        if not result.df_image.empty and "image_id" in result.df_image.columns:
            img_row = result.df_image[result.df_image["image_id"].astype(str) == str(image_id)]
            image_qc_cols = [c for c in img_row.columns if c.startswith("qc_")]
            if image_qc_cols and img_row[image_qc_cols].any(axis=1).any():
                flagged = True

        return "flagged" if flagged else "normal"

    def _populate_image_list():
        listbox.delete(0, tk.END)
        pairs = state.all_pairs()
        listbox._pairs = pairs
        for p in pairs:
            status = _image_qc_status(p["id"])
            listbox.insert(tk.END, f"{QC_BADGE_GLYPH[status]}  {p['id']}  [{p.get('group', '')}]")

    def _boundary_colors():
        if not qc_color_mode.get() or current["df_axons"] is None or current["df_axons"].empty:
            return None
        df = current["df_axons"]
        qc_cols = [c for c in df.columns if c.startswith("qc_")]
        colors = {}
        for _, row in df.iterrows():
            if bool(row.get("excluded_from_analysis", False)):
                status = "excluded"
            elif qc_cols and any(bool(row[c]) for c in qc_cols):
                status = "flagged"
            else:
                status = "normal"
            colors[int(row["axon_id"])] = QC_BOUNDARY_COLOR[status]
        return colors

    def _show_inspector(row):
        lines = [f"axon_id: {int(row['axon_id'])}"]
        for col in ("axon_area_um2", "fiber_area_um2", "myelin_area_um2", "g_ratio",
                    "circularity", "mito_count", "mito_occupancy_ratio"):
            if col in row.index:
                v = row[col]
                lines.append(f"{col}: {v:.4g}" if isinstance(v, float) else f"{col}: {v}")
        qc_hits = [c for c in row.index if c.startswith("qc_") and bool(row[c])]
        lines.append(f"qc flags: {', '.join(qc_hits) if qc_hits else 'none'}")
        if "excluded_from_analysis" in row.index:
            lines.append(f"excluded: {bool(row['excluded_from_analysis'])}")
            if row.get("exclusion_reason"):
                lines.append(f"reason: {row['exclusion_reason']}")
        inspector_var.set("\n".join(lines))

    def _on_canvas_click(event):
        df = current["df_axons"]
        if event.inaxes is None or df is None or df.empty or event.xdata is None:
            return
        dx = df["centroid_x_px"] - event.xdata
        dy = df["centroid_y_px"] - event.ydata
        dist = np.sqrt(dx ** 2 + dy ** 2)
        _show_inspector(df.loc[dist.idxmin()])

    def _render():
        if current["df_axons"] is None:
            return
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from temptation import plotting

        # Close the previous figure via pyplot BEFORE destroying its Tk
        # widgets, not after -- matplotlib's Tk backend manager tears
        # down its own canvas/toolbar association on close(), and doing
        # that against an already-.destroy()'d widget is exactly the
        # kind of ordering bug that manifests as "the second time I do
        # this, nothing happens" rather than a clean crash.
        if canvas_holder["fig"] is not None:
            import matplotlib.pyplot as plt
            plt.close(canvas_holder["fig"])
        for w in canvas_frame.winfo_children():
            w.destroy()

        fig = plotting.overlay_figure(
            current["tem"], current["mask"], current["labels_ws"], current["df_axons"],
            myelin_val=state.seg_cfg.myelin_val, axoplasm_val=state.seg_cfg.axoplasm_val,
            mito_val=state.seg_cfg.mito_val,
            title=str(current["pair"]["id"]),
            show_myelin=layer_vars["show_myelin"].get(),
            show_axoplasm=layer_vars["show_axoplasm"].get(),
            show_mito=layer_vars["show_mito"].get(),
            show_boundaries=layer_vars["show_boundaries"].get(),
            show_axon_ids=layer_vars["show_axon_ids"].get(),
            boundary_colors=_boundary_colors(),
            df_mito=current["df_mito"] if show_mito_sublayer.get() else None,
            figsize=(9, 5.5),
        )

        canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        toolbar_frame = tk.Frame(canvas_frame, bg=CARD_BG)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()

        canvas.mpl_connect("button_press_event", _on_canvas_click)

        canvas_holder["canvas"] = canvas
        canvas_holder["toolbar"] = toolbar
        canvas_holder["fig"] = fig

    def _load_and_render(pair):
        try:
            tem = dataio.read_image(pair["tem_path"])
            mask = dataio.read_mask(pair["mask_path"])
        except Exception as exc:
            inspector_var.set(f"Could not read image/mask: {exc}")
            return

        pixel_um = state.pixel_size_um or 1.0
        df_axons, _df_image, labels_ws, _resolved_mode, df_mito = pipeline.analyze_image_legacy(
            tem, mask, pixel_um, state.seg_cfg,
            qc_thresholds=state.qc_thresholds,
            exclude_qc_failed=state.exclude_qc_failed,
        )
        current.update(pair=pair, df_axons=df_axons, df_mito=df_mito,
                        tem=tem, mask=mask, labels_ws=labels_ws)
        inspector_var.set(f"{len(df_axons)} axon(s) in {pair['id']}. Click one to inspect.")
        _render()

    def _on_select(_event=None):
        sel = listbox.curselection()
        if not sel:
            return
        pairs = getattr(listbox, "_pairs", [])
        if sel[0] >= len(pairs):
            return
        _load_and_render(pairs[sel[0]])

    def _select_index(delta: int):
        """Mouse-only Previous/Next -- deliberately not bound to arrow
        keys via bind_all, which would hijack cursor movement in every
        Entry field elsewhere in the app. A layout-independent way to
        change images, on top of the list-visibility fix above."""
        pairs = getattr(listbox, "_pairs", [])
        if not pairs:
            return
        if current["pair"] is None:
            idx = 0
        else:
            try:
                idx = pairs.index(current["pair"])
            except ValueError:
                idx = 0
            idx = max(0, min(len(pairs) - 1, idx + delta))
        listbox.selection_clear(0, tk.END)
        listbox.selection_set(idx)
        listbox.see(idx)
        _load_and_render(pairs[idx])

    def _refresh():
        _populate_image_list()

    # ---------------------------------------------------------- wire it up
    listbox.bind("<<ListboxSelect>>", _on_select)

    ttk.Button(nav_row, text="◀ Previous", style="Ghost.TButton",
               command=lambda: _select_index(-1)).pack(side="left", padx=(0, 6))
    ttk.Button(nav_row, text="Next ▶", style="Ghost.TButton",
               command=lambda: _select_index(1)).pack(side="left")

    for text, var in [
        ("Myelin", layer_vars["show_myelin"]), ("Axoplasm", layer_vars["show_axoplasm"]),
        ("Mito", layer_vars["show_mito"]), ("Boundaries", layer_vars["show_boundaries"]),
        ("Axon IDs", layer_vars["show_axon_ids"]),
        ("QC coloring", qc_color_mode), ("Mito markers", show_mito_sublayer),
    ]:
        ttk.Checkbutton(toggles_row, text=text, variable=var, style="Card.TCheckbutton",
                         command=_render).pack(side="left", padx=(0, 10))

    state.on_change(_refresh)
    _refresh()

    # Test hook: Tk's <<ListboxSelect>> virtual event only fires from a
    # real mouse click, not from selection_set()/event_generate() in a
    # headless test -- stashing the loader lets tests drive the same
    # code path the click handler calls, the same pattern the legacy
    # GUI used for outer._tree in _make_table.
    parent._load_and_render = _load_and_render
