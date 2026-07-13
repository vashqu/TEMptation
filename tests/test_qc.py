"""Unit tests for temptation/qc.py (Phase 5).

Cases per IMPLEMENTATION_BLUEPRINT.md Sec 9 Phase 5 Checks / CLAUDE.md
Sec 14 case 4: g_ratio=0.97 -> qc_g_ratio_high and qc_extreme_g_ratio,
not qc_g_ratio_low; two simultaneous failures -> exclusion_reason ==
"g_ratio_high;low_circularity"; None threshold -> column all-False; row
count invariant to exclusion.
"""

import numpy as np
import pandas as pd
import pytest

from temptation.config import QCThresholds
from temptation.qc import apply_exclusions, compute_qc_flags, qc_report


def _axon_row(**overrides):
    row = {
        "axon_id": 1,
        "g_ratio": 0.65,
        "circularity": 0.7,
        "axon_area_um2": 1.0,
        "fiber_area_um2": 2.0,
        "myelin_area_um2": 0.5,
        "mito_count": 2,
    }
    row.update(overrides)
    return row


def test_extreme_g_ratio_high_not_low():
    """CLAUDE.md Sec 14 case 4: g_ratio=0.97 triggers qc_g_ratio_high and
    qc_extreme_g_ratio, never qc_g_ratio_low."""
    df = pd.DataFrame([_axon_row(g_ratio=0.97)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged["qc_g_ratio_high"].iloc[0] == True  # noqa: E712
    assert flagged["qc_extreme_g_ratio"].iloc[0] == True  # noqa: E712
    assert flagged["qc_g_ratio_low"].iloc[0] == False  # noqa: E712


def test_g_ratio_low_flag():
    df = pd.DataFrame([_axon_row(g_ratio=0.2)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged["qc_g_ratio_low"].iloc[0]
    assert flagged["qc_extreme_g_ratio"].iloc[0]
    assert not flagged["qc_g_ratio_high"].iloc[0]


def test_normal_g_ratio_no_flags():
    df = pd.DataFrame([_axon_row(g_ratio=0.65)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert not flagged["qc_g_ratio_low"].iloc[0]
    assert not flagged["qc_g_ratio_high"].iloc[0]
    assert not flagged["qc_extreme_g_ratio"].iloc[0]


def test_nan_g_ratio_never_flagged_extreme():
    """A pathological axon's NaN g_ratio must not trigger either bound --
    NaN comparisons are False, not an error or a spurious True."""
    df = pd.DataFrame([_axon_row(g_ratio=np.nan)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert not flagged["qc_g_ratio_low"].iloc[0]
    assert not flagged["qc_g_ratio_high"].iloc[0]
    assert not flagged["qc_extreme_g_ratio"].iloc[0]


def test_low_circularity_flag():
    df = pd.DataFrame([_axon_row(circularity=0.1)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged["qc_low_circularity"].iloc[0]


def test_invalid_area_relation_flag():
    df = pd.DataFrame([_axon_row(axon_area_um2=5.0, fiber_area_um2=3.0)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged["qc_invalid_area_relation"].iloc[0]


def test_invalid_area_relation_is_strict_not_equal():
    """Live-data finding: axon_area_um2 == fiber_area_um2 is the
    structurally guaranteed, EXPECTED case whenever zero myelin was
    captured for a fiber -- every pathological axon (fiber_area_px is
    defined to equal axon_area_px when watershed has no myelin channel),
    plus any 'normal'-resolved image with detached, unreachable myelin
    (F3/F4). This is already tracked via qc_no_myelin (deliberately
    non-excludable); qc_invalid_area_relation must not re-flag the same
    situation with a flag that IS excludable, or --exclude-qc-failed
    silently wipes out every axon in affected images (verified live on
    pathological_data/mask_path_1.tif before this fix -- 7/7 axons
    excluded via this flag alone)."""
    df = pd.DataFrame([_axon_row(axon_area_um2=2.0, fiber_area_um2=2.0)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert not flagged["qc_invalid_area_relation"].iloc[0]


def test_invalid_area_relation_still_fires_when_axon_exceeds_fiber():
    """The genuine anomaly this flag exists to catch (only reachable in
    practice under --assign-detached-myelin nearest) must still fire."""
    df = pd.DataFrame([_axon_row(axon_area_um2=5.0, fiber_area_um2=3.0)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged["qc_invalid_area_relation"].iloc[0]


def test_none_threshold_flag_is_always_false():
    """CLAUDE.md Sec 7: a None threshold means 'not calibrated', not
    'zero' -- the flag must be all-False, never computed against 0."""
    df = pd.DataFrame([
        _axon_row(myelin_area_um2=0.0),
        _axon_row(myelin_area_um2=-5.0),  # even a nonsensical negative value
    ])
    thresholds = QCThresholds(min_myelin_area_um2=None)
    flagged = compute_qc_flags(df, thresholds)
    assert not flagged["qc_low_myelin_area"].any()


def test_configured_threshold_flag_fires():
    df = pd.DataFrame([_axon_row(myelin_area_um2=0.05)])
    thresholds = QCThresholds(min_myelin_area_um2=0.1)
    flagged = compute_qc_flags(df, thresholds)
    assert flagged["qc_low_myelin_area"].iloc[0]


def test_no_mitochondria_and_low_count_clustering_flags():
    df = pd.DataFrame([_axon_row(mito_count=0), _axon_row(mito_count=2), _axon_row(mito_count=5)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert list(flagged["qc_no_mitochondria"]) == [True, False, False]
    assert list(flagged["qc_low_mito_count_clustering"]) == [True, True, False]


def test_no_myelin_flag_from_nan():
    df = pd.DataFrame([_axon_row(myelin_area_um2=np.nan), _axon_row(myelin_area_um2=0.5)])
    flagged = compute_qc_flags(df, QCThresholds())
    assert list(flagged["qc_no_myelin"]) == [True, False]


def test_compute_qc_flags_empty_input_no_crash():
    df = pd.DataFrame(columns=["axon_id", "g_ratio", "circularity", "axon_area_um2",
                                "fiber_area_um2", "myelin_area_um2", "mito_count"])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged.empty
    from temptation.qc import QC_FLAG_COLUMNS
    for col in QC_FLAG_COLUMNS:
        assert col in flagged.columns


def test_compute_qc_flags_never_mutates_input():
    df = pd.DataFrame([_axon_row()])
    original_cols = list(df.columns)
    compute_qc_flags(df, QCThresholds())
    assert list(df.columns) == original_cols


# ---------------------------------------------------------------------- exclusion

def test_exclusion_disabled_marks_nothing():
    df = pd.DataFrame([_axon_row(g_ratio=0.97)])
    flagged = compute_qc_flags(df, QCThresholds())
    result = apply_exclusions(flagged, enabled=False)
    assert result["excluded_from_analysis"].iloc[0] == False  # noqa: E712
    assert result["exclusion_reason"].iloc[0] == ""


def test_exclusion_enabled_single_reason():
    df = pd.DataFrame([_axon_row(g_ratio=0.97)])
    flagged = compute_qc_flags(df, QCThresholds())
    result = apply_exclusions(flagged, enabled=True)
    assert result["excluded_from_analysis"].iloc[0] == True  # noqa: E712
    assert result["exclusion_reason"].iloc[0] == "extreme_g_ratio;g_ratio_high"


def test_exclusion_reason_sorted_semicolon_joined_two_failures():
    """Blueprint's exact example: g_ratio_high;low_circularity."""
    df = pd.DataFrame([_axon_row(g_ratio=0.97, circularity=0.1)])
    flagged = compute_qc_flags(df, QCThresholds())
    result = apply_exclusions(flagged, enabled=True)
    reason = result["exclusion_reason"].iloc[0]
    tokens = reason.split(";")
    assert tokens == sorted(tokens)
    assert "g_ratio_high" in tokens
    assert "low_circularity" in tokens


def test_exclusion_never_drops_rows():
    """Row count must be invariant to --exclude-qc-failed."""
    df = pd.DataFrame([_axon_row(g_ratio=0.97), _axon_row(g_ratio=0.65)])
    flagged = compute_qc_flags(df, QCThresholds())
    result_off = apply_exclusions(flagged, enabled=False)
    result_on = apply_exclusions(flagged, enabled=True)
    assert len(result_off) == len(df)
    assert len(result_on) == len(df)


def test_qc_no_myelin_never_causes_exclusion():
    """The Sec 6.5 trap: an axon with only qc_no_myelin/qc_no_mitochondria/
    qc_low_mito_count_clustering flagged (all non-excludable) must NOT be
    excluded even when exclusion is enabled."""
    df = pd.DataFrame([_axon_row(
        g_ratio=np.nan, circularity=0.7, myelin_area_um2=np.nan, mito_count=0,
    )])
    flagged = compute_qc_flags(df, QCThresholds())
    assert flagged["qc_no_myelin"].iloc[0]
    assert flagged["qc_no_mitochondria"].iloc[0]
    assert flagged["qc_low_mito_count_clustering"].iloc[0]
    result = apply_exclusions(flagged, enabled=True)
    assert result["excluded_from_analysis"].iloc[0] == False  # noqa: E712
    assert result["exclusion_reason"].iloc[0] == ""


def test_apply_exclusions_empty_input_no_crash():
    df = pd.DataFrame(columns=["axon_id"])
    result = apply_exclusions(df, enabled=True)
    assert result.empty
    assert "excluded_from_analysis" in result.columns
    assert "exclusion_reason" in result.columns


def test_qc_report_slim_view():
    df = pd.DataFrame([{
        "image_id": "163", "group": "normal", "axon_id": 1,
        "qc_g_ratio_low": False, "qc_g_ratio_high": True,
        "excluded_from_analysis": True, "exclusion_reason": "g_ratio_high",
        "unrelated_metric_col": 42.0,
    }])
    report = qc_report(df)
    assert "unrelated_metric_col" not in report.columns
    assert "qc_g_ratio_high" in report.columns
    assert "image_id" in report.columns
