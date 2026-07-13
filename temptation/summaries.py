"""Image-level and group-level aggregation, computed over axon-level rows.

Phase 4 implements the distribution-statistics block (CLAUDE.md Sec 5.1/
5.6) and image_mvf/demyelination_index. Functions here take an explicit
set of axon rows to aggregate over -- Phase 5 will pass a QC-filtered
subset; until then, callers pass the full per-image df_axons (equivalent
to "all axons are valid"), matching IMPLEMENTATION_BLUEPRINT.md Sec 9
Phase 4's instruction to build this plumbing ahead of the QC phase rather
than retrofit it later.
"""

import numpy as np
import pandas as pd

from .mathutils import distribution_stats, safe_divide

# CLAUDE.md Sec 5.1 explicitly requests g_ratio, axon_area_um2,
# mito_density_per_um2; IMPLEMENTATION_BLUEPRINT.md Sec 4.3 extends this
# to mito_occupancy_ratio.
DISTRIBUTION_VARIABLES = ("g_ratio", "axon_area_um2", "mito_density_per_um2", "mito_occupancy_ratio")


def image_distribution_columns(df_axons_valid: pd.DataFrame, variables=DISTRIBUTION_VARIABLES) -> dict:
    """mean/std/cv/median/min/max/iqr for each variable, over the given
    axon rows, keyed as image_{stat}_{variable} (CLAUDE.md Sec 5.1/5.6
    naming -- note this is stat-before-variable, the opposite order from
    mathutils.distribution_stats's own {variable}_{stat} keys, so this
    function remaps rather than calling distribution_stats with the final
    column name as its prefix).

    Missing column or empty input -> every stat NaN for that variable
    (distribution_stats already guarantees this; never raises).
    """
    out = {}
    for var in variables:
        values = df_axons_valid[var].to_numpy() if var in df_axons_valid.columns else np.array([])
        stats = distribution_stats(values, "v")
        for tmp_key, val in stats.items():
            stat_name = tmp_key[len("v_"):]
            out[f"image_{stat_name}_{var}"] = val
    return out


def image_mvf(df_axons_valid: pd.DataFrame) -> float:
    """Myelin volume fraction over the measured fibers in this image:
    sum(myelin_area_um2) / sum(fiber_area_um2). NaN if there is no valid
    fiber area, or if every axon's myelin_area_um2 is NaN (e.g. an
    all-pathological image, where myelin is always unmeasurable -- F4)."""
    if df_axons_valid.empty or "myelin_area_um2" not in df_axons_valid.columns:
        return np.nan
    myelin_sum = df_axons_valid["myelin_area_um2"].sum(skipna=True)
    fiber_sum = df_axons_valid["fiber_area_um2"].sum(skipna=True)
    if df_axons_valid["myelin_area_um2"].isna().all():
        return np.nan
    return safe_divide(myelin_sum, fiber_sum)


def demyelination_index(myelin_area_fraction_of_fov: float, reference_myelin_fraction) -> float:
    """1 - myelin_area_fraction_of_fov / reference_myelin_fraction, clipped
    to [0, 1].

    NaN whenever no reference was supplied. IMPLEMENTATION_BLUEPRINT.md
    Sec 11 A4 flagged the reference value as a biological calibration
    choice this tool cannot derive from the code or data alone (both
    CLAUDE.md Sec 5.6 formulas collapse to axon_area/fiber_area, which is
    identically 1.0 for every pathological axon since myelin is detached
    -- F4). The user was asked and chose no default: the column is only
    computed when --demyelination-reference is explicitly supplied.
    """
    if reference_myelin_fraction is None:
        return np.nan
    ratio = safe_divide(myelin_area_fraction_of_fov, reference_myelin_fraction)
    if np.isnan(ratio):
        return np.nan
    return float(np.clip(1.0 - ratio, 0.0, 1.0))


def summarize_axon_aggregates(df_axons_subset: pd.DataFrame, resolved_mode: str) -> dict:
    """The axon-row-dependent portion of the image-level summary: legacy
    mean_* columns plus the Phase 4 distribution/image_mvf block. Does
    NOT include mask-derived quantities (mito_outside_frac, myelin
    totals, n_mito_assigned/unassigned) -- those don't depend on which
    axons are 'valid' and are computed once, not per subset.

    Phase 5 calls this twice per image: once on the full df_axons (the
    "_all" summary, always present) and once on the QC-valid subset (the
    main summary -- identical to "_all" when nothing is excluded, which
    is why calling this on the unfiltered df_axons must reproduce the
    exact legacy mean_* values; see tests/test_regression_golden.py).
    """
    if df_axons_subset.empty:
        summary = {
            "mean_axon_area_um2": np.nan,
            "mean_circularity": np.nan,
            "mean_convexity": np.nan,
            "mean_avf": np.nan,
            "mean_g_ratio": np.nan,
            "mean_myelin_thickness_um": np.nan,
            "mean_mvf": np.nan,
        }
    else:
        summary = {
            "mean_axon_area_um2": df_axons_subset["axon_area_um2"].mean(),
            "mean_circularity": df_axons_subset["circularity"].mean(),
            "mean_convexity": df_axons_subset["convexity"].mean(),
            "mean_avf": df_axons_subset["axon_vol_fraction"].mean(),
        }
        if resolved_mode == "normal":
            summary["mean_g_ratio"] = df_axons_subset["g_ratio"].mean()
            summary["mean_myelin_thickness_um"] = df_axons_subset["myelin_thickness_um"].mean()
            summary["mean_mvf"] = df_axons_subset["myelin_vol_fraction"].mean()
        else:
            summary["mean_g_ratio"] = np.nan
            summary["mean_myelin_thickness_um"] = np.nan
            summary["mean_mvf"] = np.nan

    summary.update(image_distribution_columns(df_axons_subset))
    summary["image_mvf"] = image_mvf(df_axons_subset)
    return summary


def summarize_groups(df_images: pd.DataFrame) -> pd.DataFrame:
    """One row per group: n_images, plus mean/std (ddof=1) of the key
    image-level metrics. Optional output -- the --write-group-csv CLI
    flag and group_metrics.csv file are wired up in Phase 6 ("CSV/export
    redesign"); this function is the underlying capability, usable
    standalone or from tests today.
    """
    if df_images.empty or "group" not in df_images.columns:
        return pd.DataFrame()

    metric_cols = [
        "mean_axon_area_um2", "mean_circularity", "mean_avf", "mean_g_ratio", "mean_mvf",
        "mito_outside_frac",
        "image_mean_g_ratio", "image_cv_g_ratio", "image_median_g_ratio",
        "image_mean_axon_area_um2", "image_cv_axon_area_um2",
        "image_mean_mito_density_per_um2", "image_cv_mito_density_per_um2",
        "image_mean_mito_occupancy_ratio", "image_cv_mito_occupancy_ratio",
        "image_mvf", "image_demyelination_index",
    ]
    present = [c for c in metric_cols if c in df_images.columns]

    rows = []
    for group_label, sub in df_images.groupby("group", dropna=False):
        row = {"group": group_label, "n_images": len(sub)}
        for col in present:
            row[f"{col}_mean"] = sub[col].mean(skipna=True)
            row[f"{col}_std"] = sub[col].std(skipna=True, ddof=1)
        rows.append(row)
    return pd.DataFrame(rows)
