"""Legacy-signature shim over temptation.pipeline.

Keeps `measure_image(tem, mask, pixel_length_um, ...)` with its original
keyword arguments and its original (df_axons, df_image, labels_ws,
resolved_mode) return, so the unmodified measure_nerve_gui.py (which
imports this name from measure_nerve) keeps working through Phases 1-6.
Phase 7 rewrites the GUI to call temptation.pipeline directly and this
module is removed.

`mito_hole_handling` defaults to "fill" (the Phase 2b fix for F1) and
`mito_assignment` defaults to "centroid" (the Phase 3a fix for F8) as of
this phase -- deliberately, not just for the CLI. The GUI calls this same
function with no equivalent widgets yet (those land in Phase 7's metric-
selection panel), so leaving either default at its old value here would
silently diverge CLI and GUI defaults. Pass mito_hole_handling="legacy"
and/or mito_assignment="legacy" explicitly to reproduce pre-fix numbers.

pipeline.analyze_image_legacy() itself returns a 5-tuple as of Phase 3b
(adds a per-mitochondrion df_mito). This wrapper drops that 5th element
so the GUI's exact 4-tuple contract never changes; the CLI calls
pipeline.analyze_image_legacy() directly when it needs df_mito (see
cli.py's --write-mito-csv).
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
    mito_assignment: str = "centroid",
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
        mito_assignment=mito_assignment,
    )
    df_axons, df_image, labels_ws, resolved_mode, _df_mito = analyze_image_legacy(
        tem, mask, pixel_length_um, seg_cfg,
    )
    return df_axons, df_image, labels_ws, resolved_mode
