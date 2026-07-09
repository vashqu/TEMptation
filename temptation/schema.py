"""The CSV column contract, so output shape never depends on how many
axons a particular image happened to contain."""

import numpy as np
import pandas as pd

SCHEMA_VERSION = "1.0"

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


def conform(df: pd.DataFrame, columns) -> pd.DataFrame:
    """Reindex df to exactly `columns`, filling any missing column with NaN."""
    return df.reindex(columns=list(columns))
