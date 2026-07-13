"""Quality-control flags and exclusion marking (CLAUDE.md Sec 6-8).

Two separate concerns, two separate steps:

- FLAGGED (compute_qc_flags): every qc_* boolean column is always
  computed for every axon. A flag is an observation, not an action.
- EXCLUDED (apply_exclusions): only runs when the caller opts in
  (--exclude-qc-failed). Marks excluded_from_analysis/exclusion_reason
  but never drops a row from axons.csv -- this preserves auditability
  (CLAUDE.md Sec 8) and is what makes it possible to compute both raw
  and QC-filtered image summaries from the same axon table.

CLAUDE.md Sec 6.5's trap: qc_no_myelin, qc_no_mitochondria, and
qc_low_mito_count_clustering are flags but NEVER exclusion candidates.
On this dataset, qc_no_myelin fires on every single pathological axon
(F4 -- myelin is detached, not a data-quality problem) and
qc_low_mito_count_clustering fires whenever fewer than 3 mitochondria are
present, which is the common case. Making either of these excludable
would silently delete the pathology this tool exists to characterize
-- see EXCLUDABLE_FLAGS below, which deliberately omits them.
"""

import numpy as np
import pandas as pd

from .config import QCThresholds

# Every flag this module computes, in a fixed order (also the qc_report.csv
# column order).
QC_FLAG_COLUMNS = (
    "qc_g_ratio_low",
    "qc_g_ratio_high",
    "qc_extreme_g_ratio",
    "qc_low_circularity",
    "qc_invalid_area_relation",
    "qc_low_myelin_area",
    "qc_low_axon_area",
    "qc_low_fiber_area",
    "qc_no_mitochondria",
    "qc_low_mito_count_clustering",
    "qc_no_myelin",
)

# Flags that MAY exclude an axon when --exclude-qc-failed is set.
# qc_no_mitochondria, qc_low_mito_count_clustering, and qc_no_myelin are
# deliberately excluded from this list -- see module docstring.
EXCLUDABLE_FLAGS = (
    "qc_g_ratio_low",
    "qc_g_ratio_high",
    "qc_extreme_g_ratio",
    "qc_low_circularity",
    "qc_invalid_area_relation",
    "qc_low_myelin_area",
    "qc_low_axon_area",
    "qc_low_fiber_area",
)

# Short tokens for exclusion_reason (CLAUDE.md Sec 8 example:
# "g_ratio_high;low_circularity").
_REASON_TOKEN = {
    "qc_g_ratio_low": "g_ratio_low",
    "qc_g_ratio_high": "g_ratio_high",
    "qc_extreme_g_ratio": "extreme_g_ratio",
    "qc_low_circularity": "low_circularity",
    "qc_invalid_area_relation": "invalid_area_relation",
    "qc_low_myelin_area": "low_myelin_area",
    "qc_low_axon_area": "low_axon_area",
    "qc_low_fiber_area": "low_fiber_area",
}


def _threshold_flag(series: pd.Series, threshold, empty_index) -> pd.Series:
    """series <= threshold, or all-False if threshold is None (a None
    threshold means "not calibrated for this dataset yet" -- CLAUDE.md
    Sec 7 -- not "threshold is zero")."""
    if threshold is None:
        return pd.Series(False, index=empty_index)
    return series <= threshold


def compute_qc_flags(df_axons: pd.DataFrame, thresholds: QCThresholds) -> pd.DataFrame:
    """Returns a COPY of df_axons with QC_FLAG_COLUMNS added. Never
    mutates the input, never drops rows, never raises on empty input."""
    df = df_axons.copy()

    if df.empty:
        for col in QC_FLAG_COLUMNS:
            df[col] = pd.Series(dtype=bool)
        return df

    idx = df.index
    g = df["g_ratio"]

    df["qc_g_ratio_low"] = g < thresholds.g_ratio_min
    df["qc_g_ratio_high"] = g > thresholds.g_ratio_max
    df["qc_extreme_g_ratio"] = df["qc_g_ratio_low"] | df["qc_g_ratio_high"]
    df["qc_low_circularity"] = df["circularity"] < thresholds.axon_circularity_min

    # Strict >, not >= (CLAUDE.md Sec 7 states ">="). axon_area_um2 <=
    # fiber_area_um2 is a mathematical invariant under the default
    # assign_detached_myelin="none" (axon_mask_local is constructed as a
    # subset of the fiber's watershed region -- see pipeline.py), so
    # equality is structurally guaranteed, not anomalous. It occurs
    # whenever zero myelin was captured for a fiber -- including every
    # pathological axon (fiber_area_px is defined to equal axon_area_px
    # when watershed has no myelin channel to expand into) and any
    # "normal"-resolved image whose myelin is present but detached (F3/F4:
    # unreachable by watershed, so the same zero-myelin collapse occurs
    # even though mode != "pathological"). That case is already correctly
    # tracked and deliberately non-excludable via qc_no_myelin; using >=
    # here made this flag redundant with qc_no_myelin in exactly the one
    # case that must never exclude, except this flag (unlike qc_no_myelin)
    # WAS excludable -- verified live: it silently excluded 100% of a
    # pathological image's axons under --exclude-qc-failed. Strict >
    # keeps this flag meaningful as a genuine geometric-impossibility
    # check (axon larger than its own fiber), which can only actually
    # occur under --assign-detached-myelin nearest, where fiber_area_px
    # is computed independently rather than derived as an axon superset.
    df["qc_invalid_area_relation"] = df["axon_area_um2"] > df["fiber_area_um2"]

    df["qc_low_myelin_area"] = _threshold_flag(df["myelin_area_um2"], thresholds.min_myelin_area_um2, idx)
    df["qc_low_axon_area"] = _threshold_flag(df["axon_area_um2"], thresholds.min_axon_area_um2, idx)
    df["qc_low_fiber_area"] = _threshold_flag(df["fiber_area_um2"], thresholds.min_fiber_area_um2, idx)

    df["qc_no_mitochondria"] = df["mito_count"] == 0
    df["qc_low_mito_count_clustering"] = df["mito_count"] < thresholds.min_mito_count_for_clustering
    df["qc_no_myelin"] = df["myelin_area_um2"].isna()

    # NaN-safe comparisons: pandas/numpy already return False (not NaN)
    # for e.g. `nan < 0.3`, so a pathological axon's NaN g_ratio correctly
    # yields qc_g_ratio_low=False, qc_g_ratio_high=False rather than
    # raising or propagating NaN into a boolean column.
    for col in QC_FLAG_COLUMNS:
        df[col] = df[col].fillna(False).astype(bool)

    return df


def apply_exclusions(df_axons_flagged: pd.DataFrame, enabled: bool) -> pd.DataFrame:
    """Adds excluded_from_analysis (bool) and exclusion_reason (str,
    ';'-joined sorted tokens) to a COPY of df_axons_flagged.

    When enabled=False (the default), every row gets
    excluded_from_analysis=False and exclusion_reason="" -- the columns
    are always present (schema stability), they just mark nothing.
    Never drops a row either way (CLAUDE.md Sec 8).
    """
    df = df_axons_flagged.copy()

    if df.empty:
        df["excluded_from_analysis"] = pd.Series(dtype=bool)
        df["exclusion_reason"] = pd.Series(dtype=object)
        return df

    if not enabled:
        df["excluded_from_analysis"] = False
        df["exclusion_reason"] = ""
        return df

    present_flags = [c for c in EXCLUDABLE_FLAGS if c in df.columns]

    def _reasons_for_row(row) -> str:
        tokens = sorted(_REASON_TOKEN[c] for c in present_flags if bool(row[c]))
        return ";".join(tokens)

    reasons = df.apply(_reasons_for_row, axis=1)
    df["exclusion_reason"] = reasons
    df["excluded_from_analysis"] = reasons.astype(bool)  # non-empty string -> True
    return df


def qc_report(df_axons: pd.DataFrame) -> pd.DataFrame:
    """One row per axon: identifiers + every qc_* flag + exclusion
    columns. A slim, QC-focused view of axons.csv for quick review
    without the full metric set."""
    id_cols = [c for c in ("image_id", "group", "axon_id") if c in df_axons.columns]
    flag_cols = [c for c in QC_FLAG_COLUMNS if c in df_axons.columns]
    extra_cols = [c for c in ("excluded_from_analysis", "exclusion_reason") if c in df_axons.columns]
    return df_axons[id_cols + flag_cols + extra_cols].copy()


def exclusion_reason_counts(df_axons: pd.DataFrame) -> dict:
    """How many axons were excluded for each individual reason token,
    across the whole batch. An axon excluded for two simultaneous reasons
    contributes to both counts. Used in run_manifest.json (Phase 6) so an
    exclusion-heavy run is auditable at a glance rather than requiring a
    manual scan of exclusion_reason strings."""
    if df_axons.empty or "exclusion_reason" not in df_axons.columns:
        return {}
    counts: dict = {}
    for reason in df_axons["exclusion_reason"].dropna():
        if not reason:
            continue
        for token in reason.split(";"):
            counts[token] = counts.get(token, 0) + 1
    return counts
