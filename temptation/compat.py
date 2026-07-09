"""Legacy-signature shim over temptation.pipeline.

Keeps `measure_image(tem, mask, pixel_length_um, ...)` with its original
keyword arguments and its original (df_axons, df_image, labels_ws,
resolved_mode) return, so the unmodified measure_nerve_gui.py (which
imports this name from measure_nerve) keeps working through Phases 1-6.
Phase 7 rewrites the GUI to call temptation.pipeline directly and this
module is removed.

`mito_hole_handling` defaults to "fill" (the Phase 2b fix for F1) as of
this phase -- deliberately, not just for the CLI. The GUI calls this same
function with no --mito-hole-handling-equivalent widget yet (that lands in
Phase 7's metric-selection panel), so leaving the default at "legacy"
here would silently diverge CLI and GUI defaults. Pass
mito_hole_handling="legacy" explicitly to reproduce pre-Phase-2b numbers.
"""

import numpy as np

from .config import SegmentationConfig
from .pipeline import analyze_image_legacy


def measure_image(
    tem: np.ndarray,
    mask: np.ndarray,
    pixel_length_um: float,
    myelin_val: int = 64,
    axoplasm_val: int = 192,
    mito_val: int = 128,
    mode: str = "auto",
    myelin_threshold_px: int = 200,
    smoothing_radius_px: int = 1,
    min_axon_area_px: int = 200,
    min_myelin_area_px: int = 300,
    watershed_mode: str = "weighted",
    watershed_weight: str = "radius",
    watershed_compactness: float = 0.001,
    watershed_beta: float = 1.0,
    assign_detached_myelin: str = "none",
    mito_hole_handling: str = "fill",
):
    seg_cfg = SegmentationConfig(
        myelin_val=myelin_val,
        axoplasm_val=axoplasm_val,
        mito_val=mito_val,
        mode=mode,
        myelin_threshold_px=myelin_threshold_px,
        smoothing_radius_px=smoothing_radius_px,
        min_axon_area_px=min_axon_area_px,
        min_myelin_area_px=min_myelin_area_px,
        watershed_mode=watershed_mode,
        watershed_weight=watershed_weight,
        watershed_compactness=watershed_compactness,
        watershed_beta=watershed_beta,
        assign_detached_myelin=assign_detached_myelin,
        mito_hole_handling=mito_hole_handling,
    )
    return analyze_image_legacy(tem, mask, pixel_length_um, seg_cfg)
