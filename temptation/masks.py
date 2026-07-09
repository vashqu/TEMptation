"""Turn a raw label mask into clean, mutually exclusive boolean compartments.

This is where defect F1 (IMPLEMENTATION_BLUEPRINT.md Sec 0) lives.
`mito_hole_handling="legacy"` reproduces the mitochondrial-rim bug exactly
(dilating axon_only by 1px only reaches the outer rim of each mitochondrion
hole) and exists solely so the pre-fix numbers stay reproducible for
regression testing -- do not change that branch.

`mito_hole_handling="fill"` (Phase 2b, default from this phase onward) is
the actual fix, and it is *not* simply "fill the topological holes in
axon_only" -- that first-attempt approach was tried and measured, and
rejected. `ndi.binary_fill_holes(axon_only)` only fills a hole that is
*fully* enclosed by axon_only with zero leak path to any other tissue
type; scipy treats myelin, background, and mitochondria identically as
"not axon_only" for this purpose. Measured on a real image: a 41,551px
mitochondrion that plainly abuts real axon_only territory (confirmed by
direct adjacency) also touches myelin at its far edge -- entirely normal,
since mitochondria are not restricted to the geometric center of an
axon -- and that single myelin contact was enough to disqualify the whole
component from being "filled" at all. Across that image, fill-holes
recovered only 16.5% of true mitochondrial area, barely better than the
legacy bug's 2.2%.

The correct test is connectivity, not enclosure: a mitochondrion belongs
to an axon if its connected component, unioned with axon_only, touches
real axon_only pixels -- regardless of what else it also touches. This is
implemented via `_mito_in_axon_by_connectivity` below and was verified to
recover the full mitochondrion in every case checked (all 3 components in
the same test image directly touch axon_only, confirmed via ring-adjacency
before deciding on this algorithm).
"""

from dataclasses import dataclass

import numpy as np
from skimage.measure import label as sk_label
from skimage.morphology import binary_opening, binary_closing, binary_dilation, disk

from .config import SegmentationConfig


def _mito_in_axon_by_connectivity(axon_only: np.ndarray, mito: np.ndarray) -> np.ndarray:
    """A mitochondrion belongs to the axon if its connected component
    (unioned with axon_only) contains at least one axon_only pixel --
    i.e. it directly touches real axoplasm somewhere along its boundary.
    Mitochondria with zero axon_only contact anywhere (segmentation noise
    floating in myelin/background, unconnected to any axon) are excluded."""
    union = axon_only | mito
    union_lab = sk_label(union, connectivity=2)

    # Union-labels that contain at least one axon_only pixel are "in an axon".
    labels_touching_axon = set(np.unique(union_lab[axon_only]))
    labels_touching_axon.discard(0)

    in_axon_union = np.isin(union_lab, list(labels_touching_axon))
    return mito & in_axon_union


@dataclass
class Compartments:
    axon_only: np.ndarray     # raw mask == axoplasm_val, before mito-hole handling
    axoplasm: np.ndarray      # axon_only | mito_in_axon, after smoothing
    myelin: np.ndarray        # after smoothing, mutually exclusive with axoplasm
    mito: np.ndarray          # raw mask == mito_val
    mito_in_axon: np.ndarray  # mito pixels attributed to axoplasm


def build_compartments(mask: np.ndarray, seg_cfg: SegmentationConfig) -> Compartments:
    axon_only = mask == seg_cfg.axoplasm_val
    mito = mask == seg_cfg.mito_val
    myelin = mask == seg_cfg.myelin_val

    if seg_cfg.mito_hole_handling == "legacy":
        # Verbatim reproduction of measure_nerve.py:175-177 -- reaches only
        # 1px into each mitochondrion-shaped hole (F1). Do not "fix" this
        # branch; it exists specifically to stay bug-for-bug identical.
        axon_dil = binary_dilation(axon_only, disk(1))
        mito_in_axon = mito & axon_dil
    elif seg_cfg.mito_hole_handling == "fill":
        mito_in_axon = _mito_in_axon_by_connectivity(axon_only, mito)
    else:
        raise ValueError(f"Unknown mito_hole_handling: {seg_cfg.mito_hole_handling!r}")

    axoplasm = axon_only | mito_in_axon

    if seg_cfg.smoothing_radius_px > 0:
        axoplasm = binary_opening(axoplasm, disk(seg_cfg.smoothing_radius_px))
        myelin = binary_closing(myelin, disk(seg_cfg.smoothing_radius_px))

    myelin = myelin & (~axoplasm)  # ensure mutual exclusion

    return Compartments(
        axon_only=axon_only,
        axoplasm=axoplasm,
        myelin=myelin,
        mito=mito,
        mito_in_axon=mito_in_axon,
    )


def mask_sanity(mask: np.ndarray, seg_cfg: SegmentationConfig) -> dict:
    """Fraction of pixels not in {0, myelin_val, mito_val, axoplasm_val}."""
    canonical = {0, seg_cfg.myelin_val, seg_cfg.mito_val, seg_cfg.axoplasm_val}
    unique, counts = np.unique(mask, return_counts=True)
    total = int(mask.size)
    noncanonical_px = int(sum(
        c for v, c in zip(unique, counts) if int(v) not in canonical
    ))
    return {
        "n_noncanonical_px": noncanonical_px,
        "noncanonical_frac": noncanonical_px / total if total > 0 else np.nan,
    }
