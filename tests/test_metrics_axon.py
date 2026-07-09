"""Unit tests for the Phase 2 unbiased shape metrics
(IMPLEMENTATION_BLUEPRINT.md Sec 0 F5, Sec 3 items #25-27)."""

import numpy as np
import pytest
from skimage.draw import disk as skdisk
from skimage.measure import label, regionprops

from temptation.metrics_axon import shape_metrics


def _disk_region(radius):
    img = np.zeros((2 * radius + 11, 2 * radius + 11), dtype=bool)
    rr, cc = skdisk((radius + 5, radius + 5), radius)
    img[rr, cc] = True
    lab = label(img)
    region = regionprops(lab)[0]
    return region, int(img.sum())


@pytest.mark.parametrize("radius", [10, 40, 80, 160])
def test_crofton_circularity_is_near_one_for_a_disk(radius):
    region, area_px = _disk_region(radius)
    m = shape_metrics(region, area_px, pixel_length_um=1.0)
    assert 0.95 <= m["axon_circularity_crofton"] <= 1.05


def test_legacy_circularity_is_biased_low_and_size_dependent():
    """This is F5: the legacy metric is not a bug fix target, but its
    known bias must not regress further -- pin the exact documented
    value at r=40."""
    region, area_px = _disk_region(40)
    m = shape_metrics(region, area_px, pixel_length_um=1.0)
    assert abs(m["circularity"] - 0.929) < 0.01
    # Crofton must read closer to the true value of 1.0 than legacy does.
    assert abs(m["axon_circularity_crofton"] - 1.0) < abs(m["circularity"] - 1.0)


def test_irregularity_is_exact_inverse_of_circularity():
    region, area_px = _disk_region(40)
    m = shape_metrics(region, area_px, pixel_length_um=1.0)
    assert m["axon_shape_irregularity"] == 1.0 / m["circularity"]
    assert m["axon_shape_irregularity_crofton"] == 1.0 / m["axon_circularity_crofton"]


def test_irregularity_and_circularity_never_nan_for_a_real_region():
    region, area_px = _disk_region(20)
    m = shape_metrics(region, area_px, pixel_length_um=1.0)
    for key in ("circularity", "axon_circularity_crofton", "axon_shape_irregularity", "axon_shape_irregularity_crofton"):
        assert not np.isnan(m[key]), key


def test_zero_perimeter_guard_returns_nan_not_zero_division_error():
    class FakeRegion:
        perimeter = 0.0
        perimeter_crofton = 0.0
        eccentricity = 0.0
        solidity = 1.0
        convex_image = np.array([[True]])

    m = shape_metrics(FakeRegion(), axon_area_px=100, pixel_length_um=1.0)
    assert np.isnan(m["circularity"])
    assert np.isnan(m["axon_circularity_crofton"])
    assert np.isnan(m["axon_shape_irregularity"])
    assert np.isnan(m["axon_shape_irregularity_crofton"])
