"""One overlay implementation, shared by CLI and GUI.

Phase 1 deletes the ~120-line duplicate that previously lived inline in
measure_nerve_gui.py's _on_plot (IMPLEMENTATION_BLUEPRINT.md Sec 4.3,
"temptation/plotting.py" contract). The GUI version had one feature the
CLI's make_plot() lacked: a pixel-value hover readout (format_coord). That
is folded in here as `interactive_hover` (default True) since it only
affects interactive display, never the saved PNG bytes or any CSV --
enabling it for both call sites is a strict improvement, not a behavior
change relevant to backward compatibility.
"""

import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import binary_dilation, disk

MYELIN_RGBA = np.array([0.27, 0.51, 0.71, 0.45], dtype=np.float32)
AXOPLASM_RGBA = np.array([1.00, 0.60, 0.10, 0.50], dtype=np.float32)
MITO_RGBA = np.array([0.85, 0.15, 0.15, 0.80], dtype=np.float32)


def build_overlay_rgba(mask, myelin_val=64, axoplasm_val=192, mito_val=128):
    """Tissue-type overlay array (H, W, 4), colored purely from raw mask
    values (not restricted to any watershed label set)."""
    axon_only = (mask == axoplasm_val)
    mito_raw = (mask == mito_val)
    axon_dil = binary_dilation(axon_only, disk(1))
    mito_in_axon = mito_raw & axon_dil
    axoplasm_full = axon_only | mito_in_axon

    myelin_mask = (mask == myelin_val)
    axoplasm_mask = axoplasm_full.copy()
    mito_mask = ndi.binary_fill_holes(mito_raw)
    axoplasm_mask = axoplasm_mask & ~mito_mask  # red sits on top of orange

    overlay = np.zeros((*mask.shape[:2], 4), dtype=np.float32)
    overlay[myelin_mask] = MYELIN_RGBA
    overlay[axoplasm_mask] = AXOPLASM_RGBA
    overlay[mito_mask] = MITO_RGBA
    return overlay


def _make_format_coord(mask_2d, tem_2d):
    h, w = mask_2d.shape

    def _fmt(x, y):
        ix, iy = int(round(x)), int(round(y))
        if 0 <= iy < h and 0 <= ix < w:
            mv = int(mask_2d[iy, ix])
            tv = int(tem_2d[iy, ix])
            return f"x={ix}  y={iy}  mask={mv}  tem={tv}"
        return f"x={x:.0f}  y={y:.0f}"

    return _fmt


def overlay_figure(
    tem,
    mask,
    labels_ws,
    df_axons,
    myelin_val: int = 64,
    axoplasm_val: int = 192,
    mito_val: int = 128,
    title: str = "",
    interactive_hover: bool = True,
):
    """Build the 2-panel TEM + tissue-overlay figure. Returns a Figure;
    never calls plt.show(). Verbatim geometry/coloring from the original
    make_plot()."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from skimage import measure as sk_measure

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(title, fontsize=13)

    axes[0].imshow(tem, cmap="gray", interpolation="nearest")
    axes[0].set_title("TEM image")
    axes[0].axis("off")

    axes[1].imshow(tem, cmap="gray", interpolation="nearest")

    if interactive_hover:
        mask_2d = mask if mask.ndim == 2 else mask[..., 0]
        tem_2d = tem if tem.ndim == 2 else tem[..., 0]
        fmt = _make_format_coord(mask_2d, tem_2d)
        axes[0].format_coord = fmt
        axes[1].format_coord = fmt

    overlay = build_overlay_rgba(mask, myelin_val, axoplasm_val, mito_val)
    axes[1].imshow(overlay, interpolation="nearest")

    unique_labels = np.unique(labels_ws)
    unique_labels = unique_labels[unique_labels > 0]
    for lbl in unique_labels:
        fiber_bin = (labels_ws == lbl).astype(np.uint8)
        contours = sk_measure.find_contours(fiber_bin, level=0.5)
        for contour in contours:
            axes[1].plot(
                contour[:, 1], contour[:, 0],
                color="white", linewidth=0.8, alpha=0.9,
            )

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
        Patch(facecolor=MYELIN_RGBA[:3], alpha=0.8, label="Myelin"),
        Patch(facecolor=AXOPLASM_RGBA[:3], alpha=0.8, label="Axoplasm"),
        Patch(facecolor=MITO_RGBA[:3], alpha=0.8, label="Mitochondria"),
        Line2D([0], [0], color="white", linewidth=1.2, label="Fiber boundary"),
    ]
    axes[1].legend(handles=legend_elements, loc="lower right", fontsize=7, framealpha=0.6)

    axes[1].set_title("Overlay + axon IDs")
    axes[1].axis("off")

    fig.tight_layout()
    return fig


def save_overlay(fig, save_path):
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    print(f"  [plot] saved → {save_path}")


def show_overlay_nonblocking(fig):
    """CLI --plot behavior: show briefly, then close (matches original
    make_plot() tail exactly)."""
    import matplotlib.pyplot as plt
    plt.show(block=False)
    plt.pause(0.1)
    plt.close(fig)


def show_overlay_blocking(fig):
    """GUI 'Show Overlay Plot' behavior: leave the window open."""
    import matplotlib.pyplot as plt
    plt.show(block=True)
