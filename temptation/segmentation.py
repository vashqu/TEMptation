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
