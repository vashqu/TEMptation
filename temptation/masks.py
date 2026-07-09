"""Turn a raw label mask into clean, mutually exclusive boolean compartments.

This is where defect F1 (IMPLEMENTATION_BLUEPRINT.md Sec 0) lives:
`mito_hole_handling="legacy"` reproduces the mitochondrial-rim bug exactly
(dilating axon_only by 1px only reaches the outer rim of each mitochondrion
hole). It is the only handling implemented in Phase 1. "fill" is deferred to
the gated Phase 2b.
"""

from dataclasses import dataclass

import numpy as np
from skimage.morphology import binary_opening, binary_closing, binary_dilation, disk

from .config import SegmentationConfig


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
        raise NotImplementedError(
            "mito_hole_handling='fill' lands in the gated Phase 2b "
            "(see IMPLEMENTATION_BLUEPRINT.md Sec 9, Phase 2b)."
        )
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
