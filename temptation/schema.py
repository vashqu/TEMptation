"""The CSV column contract, so output shape never depends on how many
axons a particular image happened to contain."""

import numpy as np
import pandas as pd

# SCHEMA_VERSION is the baseline/legacy value (mito_hole_handling="legacy",
# bug-for-bug identical to the pre-refactor tool). SCHEMA_VERSION_FILL
# marks output produced with the Phase 2b fix active (mito_hole_handling=
# "fill", the default from Phase 2b onward) -- see resolve_schema_version().
SCHEMA_VERSION = "1.0"
SCHEMA_VERSION_FILL = "2.0"


def resolve_schema_version(mito_hole_handling: str) -> str:
    """The schema_version column must reflect which mitochondrial-area
    algorithm actually produced a given row's numbers, not a fixed
    constant -- otherwise two runs with different --mito-hole-handling
    values would be indistinguishable in the CSV itself."""
    return SCHEMA_VERSION_FILL if mito_hole_handling == "fill" else SCHEMA_VERSION

# Exact column order of legacy axons.csv, verified against a live run of
# the pre-refactor CLI (see tests/golden/).
AXON_COLUMNS_LEGACY = (
    "image_id", "mode", "axon_id",
    "axon_area_um2", "myelin_area_um2", "fiber_area_um2",
    "d_inner_um", "d_outer_um", "g_ratio", "myelin_thickness_um",
    "perimeter_um", "eccentricity", "solidity", "circularity", "convexity",
    "axon_vol_fraction", "myelin_vol_fraction",
    "centroid_x_um", "centroid_y_um", "centroid_x_px", "centroid_y_px",
    "mito_count", "mito_area_um2", "mito_density_per_um2",
    "mito_mean_circularity", "mito_mean_form_factor", "mito_mean_feret_um",
    "mito_std_area_um2", "mito_max_area_um2", "mito_cv_area", "mito_area_skewness",
    "mito_mean_dist_centroid_um", "mito_std_dist_centroid_um",
    "nearest_neighbor_um",
)

# Exact column order of legacy image_summary.csv.
IMAGE_COLUMNS_LEGACY = (
    "image_id", "n_fibers", "fov_area_um2", "fiber_density_per_mm2",
    "mean_axon_area_um2", "mean_circularity", "mean_convexity", "mean_avf",
    "mode", "mean_g_ratio", "mean_myelin_thickness_um", "mean_mvf",
    "mito_outside_area_um2", "mito_outside_frac",
    "total_myelin_area_um2", "myelin_area_fraction_of_fov", "myelin_component_count",
)

# New in Phase 2 (unbiased shape metrics, see metrics_axon.shape_metrics
# and IMPLEMENTATION_BLUEPRINT.md Sec 0 F5). Appended after the full
# legacy block rather than interspersed, so the legacy 34 keep their exact
# values *and* positions -- IMPLEMENTATION_BLUEPRINT.md Sec 5.2 describes
# the identity block (group/image_path/mask_path/pixel_size_um/
# schema_version) as *prepended*, which would shift every legacy column's
# position; that reorder, plus the --schema {1.0,2.0} flag that makes it
# opt-in, is deliberately deferred to Phase 6 ("CSV/export redesign").
# Phase 2 stays strictly additive: nothing already exported moves.
NEW_AXON_SHAPE_COLUMNS = (
    "axon_circularity_crofton",
    "axon_shape_irregularity",
    "axon_shape_irregularity_crofton",
)

AXON_COLUMNS_V2 = AXON_COLUMNS_LEGACY + NEW_AXON_SHAPE_COLUMNS

# New in Phase 3a: diagnostics from segmentation.assign_mito_to_axons
# (global mito-to-axon assignment, fixes F8). NaN when
# mito_assignment="legacy" (no global assignment computed in that mode).
NEW_IMAGE_MITO_ASSIGNMENT_COLUMNS = ("n_mito_assigned", "n_mito_unassigned")

IMAGE_COLUMNS_V2 = IMAGE_COLUMNS_LEGACY + NEW_IMAGE_MITO_ASSIGNMENT_COLUMNS

# Identity/metadata columns, appended (not prepended -- see note above) by
# cli.py after pipeline.py returns, exactly like image_id/mode already are.
IDENTITY_COLUMNS = ("group", "image_path", "mask_path", "pixel_size_um", "schema_version")


def conform(df: pd.DataFrame, columns) -> pd.DataFrame:
    """Reindex df to exactly `columns`, filling any missing column with NaN."""
    return df.reindex(columns=list(columns))
