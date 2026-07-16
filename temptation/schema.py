"""The CSV column contract, so output shape never depends on how many
axons a particular image happened to contain."""

import numpy as np
import pandas as pd

from . import qc as _qc
from . import summaries as _summaries

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

# New in Phase 3b (metrics_mito.mito_burden_and_shape_metrics,
# metrics_spatial.mito_peripheralization_index/mito_spatial_clustering).
# Computed identically regardless of mito_assignment mode (these take a
# raw region list, not a mask -- see metrics_mito.py's module docstring),
# so unlike the Phase 3a diagnostics below, these are never NaN-by-mode.
NEW_MITO_BURDEN_SHAPE_COLUMNS = (
    "mito_total_area_um2",
    "mito_occupancy_ratio",
    "mito_mean_area_um2",
    "mito_median_area_um2",
    "mito_area_iqr",
    "mito_fragmentation_index",
    "normalized_mito_load",
    "mito_per_myelin",
    "mito_mean_aspect_ratio",
    "mito_std_aspect_ratio",
    "mito_mean_solidity",
    "mito_std_solidity",
    "mito_mean_eccentricity",
    "mito_peripheralization_index",
    "mito_mean_nn_distance_um",
    "mito_clustering_index",
)

# New in Phase 5 (qc.py): every qc_* flag is always computed
# (CLAUDE.md Sec 7), plus the exclusion-marking columns (Sec 8). Sourced
# from qc.QC_FLAG_COLUMNS directly rather than duplicated here, so the
# two lists can never drift apart.
NEW_AXON_QC_COLUMNS = _qc.QC_FLAG_COLUMNS + ("excluded_from_analysis", "exclusion_reason")

AXON_COLUMNS_V2 = (
    AXON_COLUMNS_LEGACY
    + NEW_AXON_SHAPE_COLUMNS
    + NEW_MITO_BURDEN_SHAPE_COLUMNS
    + NEW_AXON_QC_COLUMNS
)

# New in Phase 3a: diagnostics from segmentation.assign_mito_to_axons
# (global mito-to-axon assignment, fixes F8). NaN when
# mito_assignment="legacy" (no global assignment computed in that mode).
NEW_IMAGE_MITO_ASSIGNMENT_COLUMNS = ("n_mito_assigned", "n_mito_unassigned")

# New in Phase 4 (summaries.image_distribution_columns): mean/std/cv/
# median/min/max/iqr for each of summaries.DISTRIBUTION_VARIABLES, named
# image_{stat}_{variable} per CLAUDE.md Sec 5.1/5.6. Order matches
# mathutils.distribution_stats's key order (mean, std, cv, median, min,
# max, iqr).
_DISTRIBUTION_STATS_ORDER = ("mean", "std", "cv", "median", "min", "max", "iqr")
_DISTRIBUTION_VARIABLES = ("g_ratio", "axon_area_um2", "mito_density_per_um2", "mito_occupancy_ratio")
NEW_IMAGE_DISTRIBUTION_COLUMNS = tuple(
    f"image_{stat}_{var}" for var in _DISTRIBUTION_VARIABLES for stat in _DISTRIBUTION_STATS_ORDER
)

# New in Phase 4 (summaries.image_mvf, summaries.demyelination_index).
# image_demyelination_index is NaN unless --demyelination-reference is
# explicitly supplied (IMPLEMENTATION_BLUEPRINT.md Sec 11 A4 -- the
# reference is a biological calibration choice, not a default this tool
# should guess at).
NEW_IMAGE_MVF_DEMYELINATION_COLUMNS = ("image_mvf", "image_demyelination_index")

# New in Phase 5: the "_all" (unfiltered) counterpart of every key
# summaries.summarize_axon_aggregates produces, plus exclusion/flag
# diagnostics and the mask-quality check (wires up masks.mask_sanity,
# unused since Phase 2). Sourced from summarize_axon_aggregates's actual
# output keys (called once on an empty frame, which the function handles
# without raising) rather than hand-enumerated, so this can never drift
# from what pipeline.py actually emits.
_AXON_AGGREGATE_KEYS = tuple(_summaries.summarize_axon_aggregates(pd.DataFrame(), "normal").keys())
NEW_IMAGE_ALL_SUFFIXED_COLUMNS = tuple(f"{k}_all" for k in _AXON_AGGREGATE_KEYS)

NEW_IMAGE_QC_COLUMNS = (
    "image_n_axons_total",
    "image_n_axons_valid",
    "image_n_axons_excluded",
    "image_percent_flagged_axons",
    "qc_high_exclusion_rate",
    "image_noncanonical_mask_frac",
    "qc_mask_noncanonical",
)

IMAGE_COLUMNS_V2 = (
    IMAGE_COLUMNS_LEGACY
    + NEW_IMAGE_MITO_ASSIGNMENT_COLUMNS
    + NEW_IMAGE_DISTRIBUTION_COLUMNS
    + NEW_IMAGE_MVF_DEMYELINATION_COLUMNS
    + NEW_IMAGE_QC_COLUMNS
    + NEW_IMAGE_ALL_SUFFIXED_COLUMNS
)

# One row per real mitochondrion (mitochondria_metrics.csv, --write-mito-csv).
# Only produced when mito_assignment != "legacy" -- see dev/phase3a_impact.md
# for why a per-mitochondrion table requires the corrected global assignment
# to be meaningful.
MITO_COLUMNS = (
    "mito_id", "parent_axon_id",
    "mito_area_um2", "mito_perimeter_um",
    "mito_aspect_ratio", "mito_solidity", "mito_eccentricity",
    "mito_circularity", "mito_feret_um",
    "mito_centroid_x_px", "mito_centroid_y_px",
    "mito_dist_to_axon_center_um",
)

# Identity/metadata columns, appended (not prepended -- see note above) by
# cli.py after pipeline.py returns, exactly like image_id/mode already are.
IDENTITY_COLUMNS = ("group", "image_path", "mask_path", "pixel_size_um", "schema_version")


def conform(df: pd.DataFrame, columns) -> pd.DataFrame:
    """Reindex df to exactly `columns`, filling any missing column with NaN."""
    return df.reindex(columns=list(columns))
