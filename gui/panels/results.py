"""Results dashboard (blueprint Sec 8.7), the "Dashboard" tab of the
"② Review" step (alongside the per-image visual inspector in
review.py). Summary cards and a group-comparison boxplot+jitter strip,
reading already-computed columns from the last run's DataFrames;
matplotlib's own boxplot/scatter draw the statistics, nothing here
computes a metric.

The matplotlib figure is only rebuilt when state.result's identity
actually changes (tracked via `last_seen`), not on every
state.on_change() firing -- unrelated edits elsewhere in Setup (a QC
threshold keystroke, say) no longer force a matplotlib rebuild here,
since Phase-7-followup made those edits survive without nulling
`result` (see gui/state.py's mark_inputs_changed). A stale-results
banner (rather than blanking the dashboard) tells the user when the
last run's numbers may not reflect the current parameters."""

import tkinter as tk

import numpy as np

from ..widgets import CARD_BG, HDR_FG, LBL_FG, make_card

# (column, display label). Boxplot per group; suppressed (points only)
# for any group with fewer than 3 images (blueprint Sec 8.7: "a boxplot
# over two images is a lie").
AXON_LEVEL_METRICS = [
    ("g_ratio", "g-ratio"),
    ("axon_area_um2", "axon area (µm²)"),
    ("mito_occupancy_ratio", "mito occupancy"),
    ("mito_fragmentation_index", "mito fragmentation"),
]
IMAGE_LEVEL_METRICS = [
    ("image_demyelination_index", "demyelination index"),
]

MIN_IMAGES_FOR_BOX = 3


def build(parent, state):
    parent.columnconfigure(0, weight=1)
    parent.rowconfigure(2, weight=1)

    status_var = tk.StringVar(value="")
    status_label = tk.Label(parent, textvariable=status_var, bg=parent["bg"], fg=LBL_FG,
                             justify="left", anchor="w")
    status_label.grid(row=0, column=0, sticky="w", pady=(0, 10))

    cards_outer, cards_inner = make_card(parent, "Summary")
    cards_outer.grid(row=1, column=0, sticky="ew", pady=(0, 6))

    plot_outer, plot_inner = make_card(parent, "Group comparison")
    plot_outer.grid(row=2, column=0, sticky="nsew")
    plot_outer.grid_rowconfigure(0, weight=1)

    canvas_holder = {"canvas": None}
    last_seen = {"result": None}

    def _set_status(text: str):
        if text:
            status_var.set(text)
            status_label.grid()
        else:
            status_var.set("")
            status_label.grid_remove()

    def _refresh():
        result = state.result
        if result is None or result.df_axons.empty:
            _set_status("Run analysis first (① Setup) to see results here.")
            cards_outer.grid_remove()
            plot_outer.grid_remove()
            return

        cards_outer.grid()
        plot_outer.grid()
        _set_status(
            "⚠ Parameters changed since this run — showing the last completed run's results."
            if state.result_stale else ""
        )

        if result is last_seen["result"]:
            return  # only the (cheap) staleness text changed; skip the matplotlib rebuild
        last_seen["result"] = result

        for w in cards_inner.winfo_children():
            w.destroy()
        for w in plot_inner.winfo_children():
            w.destroy()
        _build_summary_cards(cards_inner, result.df_axons, result.df_image)
        _build_group_comparison(plot_inner, result.df_axons, result.df_image, canvas_holder)

    state.on_change(_refresh)
    _refresh()


def _build_summary_cards(parent, df_axons, df_image):
    n_images = df_image["image_id"].nunique() if "image_id" in df_image.columns else len(df_image)
    n_axons = len(df_axons)

    qc_cols = [c for c in df_axons.columns if c.startswith("qc_")]
    pct_flagged = 0.0
    if qc_cols and n_axons:
        pct_flagged = 100 * df_axons[qc_cols].any(axis=1).mean()

    stats = [("Images", str(n_images)), ("Axons", str(n_axons)), ("% flagged", f"{pct_flagged:.0f}%")]

    if "g_ratio" in df_axons.columns:
        g = df_axons["g_ratio"].dropna()
        if len(g):
            stats.append(("Mean g-ratio", f"{g.mean():.3f} ± {g.std():.3f}"))
    if "mito_occupancy_ratio" in df_axons.columns:
        m = df_axons["mito_occupancy_ratio"].dropna()
        if len(m):
            stats.append(("Mean mito occupancy", f"{m.mean():.3f} ± {m.std():.3f}"))

    for i, (label, value) in enumerate(stats):
        cell = tk.Frame(parent, bg=CARD_BG)
        cell.grid(row=0, column=i, sticky="w", padx=(0, 28))
        tk.Label(cell, text=value, bg=CARD_BG, fg=HDR_FG,
                 font=("TkDefaultFont", 16, "bold")).pack(anchor="w")
        tk.Label(cell, text=label, bg=CARD_BG, fg="#999999", font=("TkDefaultFont", 8)).pack(anchor="w")


def _build_group_comparison(parent, df_axons, df_image, canvas_holder):
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure

    if "group" not in df_axons.columns or df_axons["group"].dropna().empty:
        tk.Label(parent, text="No group labels present -- assign groups in ① Setup → Data.",
                 bg=CARD_BG, fg="#999999").pack(anchor="w")
        return

    groups = sorted(g for g in df_axons["group"].dropna().unique())
    if "image_id" in df_image.columns:
        n_images_per_group = df_image.groupby("group")["image_id"].nunique()
    else:
        n_images_per_group = df_image.groupby("group").size()

    metrics = [(c, lbl, df_axons) for c, lbl in AXON_LEVEL_METRICS if c in df_axons.columns]
    metrics += [(c, lbl, df_image) for c, lbl in IMAGE_LEVEL_METRICS if c in df_image.columns]

    if not metrics:
        tk.Label(parent, text="No comparable metrics in this run's output.",
                 bg=CARD_BG, fg="#999999").pack(anchor="w")
        return

    rng = np.random.default_rng(0)
    fig = Figure(figsize=(max(2.6 * len(metrics), 6), 3.4), dpi=90)
    for i, (col, label, source_df) in enumerate(metrics):
        ax = fig.add_subplot(1, len(metrics), i + 1)
        box_data, box_positions = [], []
        for gi, g in enumerate(groups):
            n_img = int(n_images_per_group.get(g, 0))
            values = source_df.loc[source_df["group"] == g, col].dropna().to_numpy()
            if len(values) == 0:
                continue
            jitter = (rng.random(len(values)) - 0.5) * 0.3
            ax.scatter(np.full(len(values), gi + 1) + jitter, values, s=8, alpha=0.5, color="#007BFF")
            if n_img >= MIN_IMAGES_FOR_BOX:
                box_data.append(values)
                box_positions.append(gi + 1)
        if box_data:
            ax.boxplot(box_data, positions=box_positions, widths=0.5, showfliers=False)
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels([f"{g}\n(n={int(n_images_per_group.get(g, 0))})" for g in groups], fontsize=7)
        ax.set_title(label, fontsize=9)
        ax.tick_params(labelsize=7)

    fig.tight_layout()

    if canvas_holder["canvas"] is not None:
        canvas_holder["canvas"].get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(fig, master=parent)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)
    canvas_holder["canvas"] = canvas
