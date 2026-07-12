"""Unit tests for the Phase 3b mitochondrial burden/shape/spatial metrics
(metrics_mito.mito_burden_and_shape_metrics,
metrics_spatial.mito_peripheralization_index,
metrics_spatial.mito_spatial_clustering).

Covers CLAUDE.md Sec 14's cases 2 (two mitochondria of known area), 3
(zero mitochondria), and 5 (clustering unstable below 3 mitochondria).
"""

import numpy as np
import pytest
from skimage.draw import disk as skdisk
from skimage.measure import label, regionprops

from temptation.metrics_mito import mito_burden_and_shape_metrics
from temptation.metrics_spatial import mito_peripheralization_index, mito_spatial_clustering


def _regions_from_disks(specs, shape=(200, 200)):
    """specs: list of (center_row, center_col, radius). Returns a list of
    regionprops, one per disk (assumes disks don't touch)."""
    img = np.zeros(shape, dtype=bool)
    for cy, cx, r in specs:
        rr, cc = skdisk((cy, cx), r, shape=shape)
        img[rr, cc] = True
    lab = label(img, connectivity=2)
    return regionprops(lab)


PX = 1.0  # keep areas in raw px^2 = um^2 for hand-checkable numbers


def test_two_mitochondria_of_known_area_count_and_fragmentation():
    """CLAUDE.md Sec 14 case 2: mito_count == 2, total == sum of the two
    known areas, fragmentation_index == 2/total."""
    regions = _regions_from_disks([(50, 50, 5), (150, 150, 8)])
    assert len(regions) == 2
    known_total_px = regions[0].area + regions[1].area

    m = mito_burden_and_shape_metrics(
        regions, axon_area_um2=1000.0, fiber_area_um2=2000.0, myelin_area_um2=500.0,
        pixel_length_um=PX,
    )
    assert m["mito_total_area_um2"] == pytest.approx(known_total_px, rel=1e-9)
    assert m["mito_fragmentation_index"] == pytest.approx(2 / known_total_px, rel=1e-9)
    assert m["mito_occupancy_ratio"] == pytest.approx(known_total_px / 1000.0, rel=1e-9)
    assert m["normalized_mito_load"] == pytest.approx(known_total_px / 2000.0, rel=1e-9)
    assert m["mito_per_myelin"] == pytest.approx(known_total_px / 500.0, rel=1e-9)


def test_zero_mitochondria_no_crash_and_correct_defaults():
    """CLAUDE.md Sec 14 case 3: no crash; density/occupancy 0; shape
    summaries NaN."""
    m = mito_burden_and_shape_metrics(
        [], axon_area_um2=1000.0, fiber_area_um2=2000.0, myelin_area_um2=500.0,
        pixel_length_um=PX,
    )
    assert m["mito_total_area_um2"] == 0.0
    assert m["mito_occupancy_ratio"] == 0.0
    assert m["normalized_mito_load"] == 0.0
    assert m["mito_per_myelin"] == 0.0
    assert np.isnan(m["mito_fragmentation_index"])
    assert np.isnan(m["mito_mean_area_um2"])
    assert np.isnan(m["mito_median_area_um2"])
    assert np.isnan(m["mito_area_iqr"])
    assert np.isnan(m["mito_mean_aspect_ratio"])
    assert np.isnan(m["mito_std_aspect_ratio"])
    assert np.isnan(m["mito_mean_solidity"])
    assert np.isnan(m["mito_std_solidity"])
    assert np.isnan(m["mito_mean_eccentricity"])


def test_zero_axon_area_gives_nan_not_zero_division_error():
    regions = _regions_from_disks([(50, 50, 5)])
    m = mito_burden_and_shape_metrics(
        regions, axon_area_um2=0.0, fiber_area_um2=0.0, myelin_area_um2=np.nan,
        pixel_length_um=PX,
    )
    assert np.isnan(m["mito_occupancy_ratio"])
    assert np.isnan(m["normalized_mito_load"])
    assert np.isnan(m["mito_per_myelin"])


def test_invalid_myelin_area_gives_nan_mito_per_myelin():
    regions = _regions_from_disks([(50, 50, 5)])
    m = mito_burden_and_shape_metrics(
        regions, axon_area_um2=1000.0, fiber_area_um2=2000.0, myelin_area_um2=np.nan,
        pixel_length_um=PX,
    )
    assert np.isnan(m["mito_per_myelin"])

    m2 = mito_burden_and_shape_metrics(
        regions, axon_area_um2=1000.0, fiber_area_um2=2000.0, myelin_area_um2=0.0,
        pixel_length_um=PX,
    )
    assert np.isnan(m2["mito_per_myelin"])


def test_single_mitochondrion_std_fields_are_nan_not_zero():
    """New (Phase 3b) metrics use ddof=1, unlike the legacy std columns
    which default to 0.0 for a single value -- single-mito std must be
    NaN here."""
    regions = _regions_from_disks([(50, 50, 5)])
    m = mito_burden_and_shape_metrics(
        regions, axon_area_um2=1000.0, fiber_area_um2=2000.0, myelin_area_um2=500.0,
        pixel_length_um=PX,
    )
    assert np.isnan(m["mito_std_aspect_ratio"])
    assert np.isnan(m["mito_std_solidity"])
    assert not np.isnan(m["mito_mean_aspect_ratio"])


def test_peripheralization_index_is_zero_at_center():
    # mean distance 0 (mitochondria exactly at axon center) -> index 0
    idx = mito_peripheralization_index(0.0, axon_area_um2=100 * np.pi)
    assert idx == 0.0


def test_peripheralization_index_scales_with_radius():
    # axon of radius 10um, mitochondria on average 10um from center -> index ~1
    axon_area_um2 = np.pi * 10 ** 2
    idx = mito_peripheralization_index(10.0, axon_area_um2)
    assert idx == pytest.approx(1.0, rel=1e-9)


def test_peripheralization_nan_for_invalid_axon_area():
    assert np.isnan(mito_peripheralization_index(5.0, 0.0))
    assert np.isnan(mito_peripheralization_index(5.0, np.nan))


def test_clustering_nan_below_min_count():
    """CLAUDE.md Sec 14 case 5: clustering_index is NaN for mito_count < 3,
    even though mito_mean_nn_distance_um (2 points) is well-defined."""
    regions = _regions_from_disks([(50, 50, 5), (150, 150, 5)])
    result = mito_spatial_clustering(regions, axon_area_um2=1000.0, pixel_length_um=PX)
    assert not np.isnan(result["mito_mean_nn_distance_um"])
    assert np.isnan(result["mito_clustering_index"])


def test_clustering_nan_for_zero_or_one_mitochondria():
    assert np.isnan(mito_spatial_clustering([], 1000.0, PX)["mito_mean_nn_distance_um"])
    one_region = _regions_from_disks([(50, 50, 5)])
    r1 = mito_spatial_clustering(one_region, 1000.0, PX)
    assert np.isnan(r1["mito_mean_nn_distance_um"])
    assert np.isnan(r1["mito_clustering_index"])


def test_clustering_index_computable_at_three_or_more():
    regions = _regions_from_disks([(30, 30, 4), (30, 170, 4), (170, 30, 4)])
    result = mito_spatial_clustering(regions, axon_area_um2=40000.0, pixel_length_um=PX)
    assert not np.isnan(result["mito_mean_nn_distance_um"])
    assert not np.isnan(result["mito_clustering_index"])
    assert result["mito_clustering_index"] > 0


def test_clustering_index_around_one_for_grid_like_spread():
    """A grid arrangement approximates CSR reasonably well over a large
    axon -- clustering index should land in a broad plausible range, not
    at an extreme."""
    specs = [(30 + 40 * i, 30 + 40 * j, 4) for i in range(4) for j in range(4)]
    regions = _regions_from_disks(specs, shape=(200, 200))
    result = mito_spatial_clustering(regions, axon_area_um2=200 * 200, pixel_length_um=PX)
    assert 0.1 < result["mito_clustering_index"] < 10
