"""Synthetic-geometry tests for the Phase 2b mitochondrial-hole fix
(IMPLEMENTATION_BLUEPRINT.md Sec 0 F1, Sec 9 Phase 2b).

These pin the exact behavior discovered while implementing the fix: a
naive `ndi.binary_fill_holes(axon_only)` approach was tried first and
rejected because a mitochondrion touching myelin (very common -- they are
not restricted to the axon's geometric interior) breaks full topological
enclosure. The shipped algorithm is connectivity-based instead.
"""

import numpy as np
from skimage.draw import disk as skdisk

from temptation.config import SegmentationConfig
from temptation.masks import build_compartments, _mito_in_axon_by_connectivity


def _blank_mask(shape=(200, 200)):
    return np.zeros(shape, dtype=np.uint8)


def _paint_disk(mask, center, radius, value):
    rr, cc = skdisk(center, radius, shape=mask.shape)
    mask[rr, cc] = value
    return mask


def test_interior_mito_hole_fully_captured():
    """A mitochondrion entirely inside an axon (myelin ring around the
    whole fiber) must be captured 100% under fill mode."""
    mask = _blank_mask()
    mask = _paint_disk(mask, (100, 100), 60, 64)   # myelin ring (fiber)
    mask = _paint_disk(mask, (100, 100), 50, 192)  # axoplasm
    mask = _paint_disk(mask, (90, 100), 8, 128)    # mito #1, interior
    mask = _paint_disk(mask, (110, 90), 5, 128)    # mito #2, interior

    true_mito_px = int((mask == 128).sum())

    cfg = SegmentationConfig(mito_hole_handling="fill")
    comps = build_compartments(mask, cfg)
    assert comps.mito_in_axon.sum() == true_mito_px

    # axoplasm area should equal axon_only + all mito (the "holes" are the
    # mito pixels themselves, not extra area beyond them)
    axon_only_px = int((mask == 192).sum())
    assert comps.axoplasm.sum() == axon_only_px + true_mito_px


def test_mito_touching_myelin_boundary_still_fully_captured():
    """This is the exact case that broke the first-attempt fill-holes
    algorithm: a mitochondrion that touches myelin at one edge (so it is
    not fully enclosed by axon_only alone) but clearly abuts real
    axoplasm elsewhere. Connectivity must still capture it in full."""
    mask = _blank_mask()
    mask = _paint_disk(mask, (100, 100), 60, 64)
    mask = _paint_disk(mask, (100, 100), 50, 192)
    # mito centered near the axon/myelin boundary (radius 50), overlapping
    # into the myelin ring
    mask = _paint_disk(mask, (100, 145), 10, 128)

    true_mito_px = int((mask == 128).sum())
    assert true_mito_px > 0

    cfg = SegmentationConfig(mito_hole_handling="fill")
    comps = build_compartments(mask, cfg)
    assert comps.mito_in_axon.sum() == true_mito_px


def test_orphaned_mito_not_touching_any_axon_is_excluded():
    """A mito-labeled blob with zero adjacency to axon_only (segmentation
    noise, or genuinely non-axonal mitochondria) must not be attributed
    to any axon."""
    mask = _blank_mask()
    mask = _paint_disk(mask, (60, 60), 30, 64)
    mask = _paint_disk(mask, (60, 60), 22, 192)
    # a mito blob far away, touching nothing
    mask = _paint_disk(mask, (170, 170), 6, 128)

    cfg = SegmentationConfig(mito_hole_handling="fill")
    comps = build_compartments(mask, cfg)
    assert comps.mito_in_axon.sum() == 0


def test_connectivity_helper_matches_build_compartments():
    mask = _blank_mask()
    mask = _paint_disk(mask, (100, 100), 60, 64)
    mask = _paint_disk(mask, (100, 100), 50, 192)
    mask = _paint_disk(mask, (90, 100), 8, 128)

    axon_only = mask == 192
    mito = mask == 128
    direct = _mito_in_axon_by_connectivity(axon_only, mito)

    cfg = SegmentationConfig(mito_hole_handling="fill")
    comps = build_compartments(mask, cfg)
    assert np.array_equal(comps.mito_in_axon, direct)


def test_legacy_mode_still_undercaptures_on_the_same_synthetic_case():
    """Confirms the legacy branch is untouched by this fix -- it must
    still only capture a thin rim, not the interior mitochondrion."""
    mask = _blank_mask()
    mask = _paint_disk(mask, (100, 100), 60, 64)
    mask = _paint_disk(mask, (100, 100), 50, 192)
    mask = _paint_disk(mask, (90, 100), 8, 128)  # radius-8 mito, interior

    true_mito_px = int((mask == 128).sum())

    cfg = SegmentationConfig(mito_hole_handling="legacy")
    comps = build_compartments(mask, cfg)
    assert 0 < comps.mito_in_axon.sum() < true_mito_px, (
        "legacy mode should still under-capture (rim only), not fully "
        "capture the mitochondrion"
    )


def test_fill_and_legacy_are_both_valid_and_produce_different_results():
    mask = _blank_mask()
    mask = _paint_disk(mask, (100, 100), 60, 64)
    mask = _paint_disk(mask, (100, 100), 50, 192)
    mask = _paint_disk(mask, (90, 100), 8, 128)

    fill_comps = build_compartments(mask, SegmentationConfig(mito_hole_handling="fill"))
    legacy_comps = build_compartments(mask, SegmentationConfig(mito_hole_handling="legacy"))
    assert fill_comps.mito_in_axon.sum() > legacy_comps.mito_in_axon.sum()
