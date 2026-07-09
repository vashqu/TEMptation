"""Mode resolution, axon labelling, and watershed. No physical units here."""

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import label, regionprops
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed

from .config import SegmentationConfig
from .masks import Compartments


def resolve_mode(mask: np.ndarray, seg_cfg: SegmentationConfig) -> str:
    """Verbatim port of measure_nerve.py:156-165.

    Uses the raw mask (not the smoothed Compartments), matching the
    original call order: mode is resolved before compartments are built.
    """
    myelin_px_count = int((mask == seg_cfg.myelin_val).sum())

    if seg_cfg.mode == "auto":
        return "normal" if myelin_px_count >= seg_cfg.myelin_threshold_px else "pathological"
    elif seg_cfg.mode in ("normal", "pathological"):
        return seg_cfg.mode
    else:
        raise ValueError(f"mode must be 'normal', 'pathological', or 'auto', got {seg_cfg.mode!r}")


def _build_weighted_distance(axon_lab: np.ndarray, weight: str, beta: float) -> np.ndarray:
    """
    Return a weighted distance map for watershed.

    elevation = max(distance - beta * weight_map, 0)

    weight='radius' : weight_map = equivalent_radius  (larger axons expand more)
    weight='area'   : weight_map = sqrt(area)          (similar idea, area-based)
    beta            : bias strength; higher = larger axons claim more territory
    """
    distance, indices = ndi.distance_transform_edt(axon_lab == 0, return_indices=True)
    nearest_axon_map = axon_lab[indices[0], indices[1]]

    props = regionprops(axon_lab)
    max_label = axon_lab.max()
    weight_lookup = np.ones(max_label + 1, dtype=np.float32)   # default 1 (safe)

    for p in props:
        if weight == "radius":
            weight_lookup[p.label] = max(p.equivalent_diameter / 2.0, 1e-6)
        elif weight == "area":
            weight_lookup[p.label] = max(np.sqrt(p.area), 1e-6)
        else:
            raise ValueError(f"Unknown watershed weight: {weight!r}")

    weight_map = weight_lookup[nearest_axon_map]
    elevation = distance - beta * weight_map
    return np.maximum(elevation, 0.0)


def label_axons(comps: Compartments, seg_cfg: SegmentationConfig) -> np.ndarray:
    """Verbatim port of measure_nerve.py:186-190."""
    axon_lab = label(comps.axoplasm, connectivity=2)
    axon_lab = remove_small_objects(axon_lab, min_size=seg_cfg.min_axon_area_px)
    axon_lab = label(axon_lab > 0, connectivity=2)
    return axon_lab


def run_watershed(
    axon_lab: np.ndarray,
    comps: Compartments,
    resolved_mode: str,
    seg_cfg: SegmentationConfig,
) -> np.ndarray:
    """Verbatim port of measure_nerve.py:192-210."""
    if resolved_mode == "normal":
        allowed_mask = (axon_lab > 0) | comps.myelin
    else:
        allowed_mask = axon_lab > 0

    if seg_cfg.watershed_mode == "weighted":
        elevation = _build_weighted_distance(axon_lab, seg_cfg.watershed_weight, beta=seg_cfg.watershed_beta)
    else:  # simple
        elevation = ndi.distance_transform_edt(axon_lab == 0)

    labels_ws = watershed(
        elevation,
        markers=axon_lab,
        mask=allowed_mask,
        compactness=seg_cfg.watershed_compactness,
    )
    return labels_ws


def nearest_axon_map(axon_lab: np.ndarray) -> np.ndarray:
    """Verbatim port of measure_nerve.py:215-219 (detached-myelin assignment map)."""
    _, nn_idx = ndi.distance_transform_edt(axon_lab == 0, return_indices=True)
    return axon_lab[nn_idx[0], nn_idx[1]]


def assign_mito_to_axons(mito_in_axon: np.ndarray, labels_ws: np.ndarray, method: str = "centroid"):
    """Label every mitochondrion in `mito_in_axon` once, globally, and
    assign each one to exactly one watershed fiber (CLAUDE.md Sec 12.4,
    IMPLEMENTATION_BLUEPRINT.md F8). This replaces the legacy per-fiber
    crop-and-intersect approach, which double-counts (and fragments the
    shape of) any mitochondrion straddling a watershed boundary between
    two axons -- verified absent in the current dataset (0/1005
    mitochondria straddle a boundary across all 99 masks) but not
    guaranteed on other data.

    method:
      "centroid" (default) -- assign to the fiber at the mitochondrion's
        centroid pixel. Preferred per CLAUDE.md Sec 12.4.
      "overlap" -- assign to the fiber with maximum pixel overlap
        (robust fallback for boundary cases where the centroid itself
        might land just outside every fiber, e.g. in myelin).

    Returns (mito_lab, assignment) where mito_lab is the global mito label
    array (connectivity=2, matching the rest of the pipeline) and
    assignment is {mito_label: fiber_label}; a mitochondrion with no
    fiber assignment (centroid lands in myelin/background and overlap is
    also zero) is omitted from `assignment`.
    """
    mito_lab = label(mito_in_axon, connectivity=2)
    assignment = {}

    for r in regionprops(mito_lab):
        if method == "centroid":
            cy, cx = r.centroid
            iy, ix = int(round(cy)), int(round(cx))
            fiber_label = int(labels_ws[iy, ix])
            if fiber_label == 0:
                # Centroid landed outside any fiber (e.g. an elongated or
                # concave mitochondrion whose centroid falls in myelin) --
                # fall back to max overlap for this component only.
                fiber_label = _max_overlap_label(mito_lab, r.label, labels_ws)
        elif method == "overlap":
            fiber_label = _max_overlap_label(mito_lab, r.label, labels_ws)
        else:
            raise ValueError(f"Unknown mito_assignment method: {method!r}")

        if fiber_label != 0:
            assignment[r.label] = fiber_label

    return mito_lab, assignment


def _max_overlap_label(mito_lab: np.ndarray, mito_label: int, labels_ws: np.ndarray) -> int:
    ys, xs = np.where(mito_lab == mito_label)
    ws_here = labels_ws[ys, xs]
    ws_here = ws_here[ws_here != 0]
    if ws_here.size == 0:
        return 0
    values, counts = np.unique(ws_here, return_counts=True)
    return int(values[np.argmax(counts)])
