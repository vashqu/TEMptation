"""Unit tests for temptation/summaries.py (Phase 4).

Hand-checkable cases per IMPLEMENTATION_BLUEPRINT.md Sec 10.2/Sec 9 Phase
4 Checks: [1,2,3,4] -> median=2.5, iqr=1.5, std(ddof=1)=1.29099...,
cv=0.51639...; empty input -> all NaN, no exception; single value ->
std=NaN, cv=NaN, iqr=0.
"""

import numpy as np
import pandas as pd
import pytest

from temptation.summaries import (
    demyelination_index,
    image_distribution_columns,
    image_mvf,
    summarize_groups,
)


def test_image_distribution_columns_hand_checked_values():
    df = pd.DataFrame({"g_ratio": [1, 2, 3, 4]})
    out = image_distribution_columns(df, variables=("g_ratio",))
    assert out["image_mean_g_ratio"] == pytest.approx(2.5)
    assert out["image_median_g_ratio"] == pytest.approx(2.5)
    assert out["image_std_g_ratio"] == pytest.approx(1.2909944487358056)
    assert out["image_cv_g_ratio"] == pytest.approx(0.5163977794943222)
    assert out["image_iqr_g_ratio"] == pytest.approx(1.5)
    assert out["image_min_g_ratio"] == 1.0
    assert out["image_max_g_ratio"] == 4.0


def test_image_distribution_columns_naming_is_stat_before_variable():
    """CLAUDE.md Sec 5.6 naming: image_{stat}_{variable}, e.g.
    image_mean_g_ratio -- not image_g_ratio_mean."""
    df = pd.DataFrame({"g_ratio": [1, 2, 3, 4]})
    out = image_distribution_columns(df, variables=("g_ratio",))
    assert "image_mean_g_ratio" in out
    assert "image_g_ratio_mean" not in out


def test_image_distribution_columns_empty_input_all_nan_no_exception():
    df = pd.DataFrame({"g_ratio": []})
    out = image_distribution_columns(df, variables=("g_ratio",))
    assert all(np.isnan(v) for v in out.values())


def test_image_distribution_columns_missing_column_all_nan():
    df = pd.DataFrame({"other_col": [1, 2, 3]})
    out = image_distribution_columns(df, variables=("g_ratio",))
    assert all(np.isnan(v) for v in out.values())


def test_image_distribution_columns_single_value():
    df = pd.DataFrame({"g_ratio": [7.0]})
    out = image_distribution_columns(df, variables=("g_ratio",))
    assert np.isnan(out["image_std_g_ratio"])
    assert np.isnan(out["image_cv_g_ratio"])
    assert out["image_iqr_g_ratio"] == 0.0
    assert out["image_mean_g_ratio"] == 7.0


def test_image_distribution_columns_multiple_variables():
    df = pd.DataFrame({"g_ratio": [1, 2, 3], "axon_area_um2": [10, 20, 30]})
    out = image_distribution_columns(df, variables=("g_ratio", "axon_area_um2"))
    assert "image_mean_g_ratio" in out
    assert "image_mean_axon_area_um2" in out
    assert out["image_mean_axon_area_um2"] == pytest.approx(20.0)


def test_image_mvf_basic():
    df = pd.DataFrame({
        "myelin_area_um2": [1.0, 2.0, 3.0],
        "fiber_area_um2": [4.0, 5.0, 6.0],
    })
    assert image_mvf(df) == pytest.approx(6.0 / 15.0)


def test_image_mvf_all_nan_myelin_returns_nan():
    """F4: an all-pathological image has myelin_area_um2 = NaN for every
    axon; image_mvf must be NaN, not a division involving 0."""
    df = pd.DataFrame({
        "myelin_area_um2": [np.nan, np.nan],
        "fiber_area_um2": [4.0, 5.0],
    })
    assert np.isnan(image_mvf(df))


def test_image_mvf_empty_dataframe():
    assert np.isnan(image_mvf(pd.DataFrame({"myelin_area_um2": [], "fiber_area_um2": []})))


def test_demyelination_index_nan_without_reference():
    """The user explicitly chose no default: NaN unless a reference is
    supplied."""
    assert np.isnan(demyelination_index(0.10, None))


def test_demyelination_index_computed_with_reference():
    # myelin fraction = 0.10, reference = 0.30 -> 1 - 0.10/0.30 = 0.6667
    idx = demyelination_index(0.10, 0.30)
    assert idx == pytest.approx(1.0 - 0.10 / 0.30)


def test_demyelination_index_clipped_to_zero_when_above_reference():
    # myelin fraction exceeds reference -> ratio > 1 -> 1-ratio < 0 -> clip to 0
    idx = demyelination_index(0.50, 0.30)
    assert idx == 0.0


def test_demyelination_index_clipped_to_one_for_near_zero_myelin():
    idx = demyelination_index(0.0001, 0.30)
    assert 0.99 < idx <= 1.0


def test_demyelination_index_nan_for_zero_reference():
    assert np.isnan(demyelination_index(0.10, 0.0))


def test_summarize_groups_basic():
    df = pd.DataFrame({
        "group": ["normal", "normal", "pathological"],
        "mean_g_ratio": [0.6, 0.7, np.nan],
        "image_mean_g_ratio": [0.6, 0.7, 0.5],
    })
    out = summarize_groups(df)
    assert set(out["group"]) == {"normal", "pathological"}
    normal_row = out[out.group == "normal"].iloc[0]
    assert normal_row["n_images"] == 2
    assert normal_row["mean_g_ratio_mean"] == pytest.approx(0.65)
    assert normal_row["image_mean_g_ratio_mean"] == pytest.approx(0.65)


def test_summarize_groups_empty_input():
    assert summarize_groups(pd.DataFrame()).empty


def test_summarize_groups_no_group_column():
    assert summarize_groups(pd.DataFrame({"x": [1, 2]})).empty
