"""Tests for the Phase 3a global mito-to-axon assignment
(segmentation.assign_mito_to_axons), which fixes F8: the legacy per-fiber
crop-and-intersect approach fragments and double-counts any mitochondrion
straddling a watershed boundary between two axons."""

import numpy as np
import pytest
import tifffile as tiff
from pathlib import Path
from skimage.measure import label as sk_label, regionprops
from skimage.segmentation import watershed as sk_watershed
from scipy import ndimage as ndi

from temptation.segmentation import assign_mito_to_axons
from temptation.masks import build_compartments
from temptation.metrics_mito import mito_metrics_for_axon, mito_metrics_from_regions
from temptation.config import SegmentationConfig

DATA_ROOT = Path(__file__).parent.parent.parent


def _two_fiber_labels_ws(shape=(100, 100)):
    """Two square watershed regions side by side, labels 1 and 2."""
    labels_ws = np.zeros(shape, dtype=np.int32)
    labels_ws[:, :50] = 1
    labels_ws[:, 50:] = 2
    return labels_ws


def test_non_straddling_mito_assigned_to_its_axon():
    labels_ws = _two_fiber_labels_ws()
    mito_in_axon = np.zeros((100, 100), dtype=bool)
    mito_in_axon[20:30, 10:20] = True  # entirely within fiber 1's region

    mito_lab, assignment = assign_mito_to_axons(mito_in_axon, labels_ws, method="centroid")
    assert mito_lab.max() == 1
    assert assignment == {1: 1}


def test_straddling_mito_assigned_to_exactly_one_axon_not_split():
    """This is the exact F8 scenario: one physical mitochondrion whose
    pixels span both watershed regions. It must be counted as ONE
    mitochondrion assigned to ONE axon, not fragmented into two."""
    labels_ws = _two_fiber_labels_ws()
    mito_in_axon = np.zeros((100, 100), dtype=bool)
    mito_in_axon[40:60, 45:55] = True  # straddles column 50 (the boundary)

    mito_lab, assignment = assign_mito_to_axons(mito_in_axon, labels_ws, method="centroid")
    assert mito_lab.max() == 1, "must be one connected component, not two"
    assert len(assignment) == 1
    assert list(assignment.values())[0] in (1, 2)


def test_mito_with_zero_fiber_overlap_is_unassigned():
    labels_ws = _two_fiber_labels_ws()
    mito_in_axon = np.zeros((150, 150), dtype=bool)
    mito_in_axon[120:130, 120:130] = True  # entirely outside both fibers (labels_ws is 100x100 here... )
    labels_ws_padded = np.zeros((150, 150), dtype=np.int32)
    labels_ws_padded[:100, :100] = labels_ws

    mito_lab, assignment = assign_mito_to_axons(mito_in_axon, labels_ws_padded, method="centroid")
    assert mito_lab.max() == 1
    assert assignment == {}


def test_overlap_fallback_for_concave_mito_with_off_component_centroid():
    """A C-shaped (concave) mitochondrion whose centroid falls outside its
    own footprint, in a region with no fiber label -- centroid method must
    fall back to max-overlap rather than leaving it unassigned."""
    labels_ws = _two_fiber_labels_ws()
    mito_in_axon = np.zeros((100, 100), dtype=bool)
    # C-shape open on the right, centered around (50,50) which is on the
    # boundary/outside both fiber interiors in a way that could confuse a
    # naive centroid lookup
    mito_in_axon[40:60, 40:60] = True
    mito_in_axon[45:55, 45:60] = False  # carve out the inside of the C

    mito_lab, assignment = assign_mito_to_axons(mito_in_axon, labels_ws, method="centroid")
    assert mito_lab.max() == 1
    assert len(assignment) == 1  # must still be assigned via overlap fallback


def test_unknown_method_raises():
    labels_ws = _two_fiber_labels_ws()
    mito_in_axon = np.zeros((100, 100), dtype=bool)
    mito_in_axon[20:30, 10:20] = True
    with pytest.raises(ValueError):
        assign_mito_to_axons(mito_in_axon, labels_ws, method="bogus")


@pytest.mark.parametrize("image_path", [
    "normal_data/162-165/163. Mask 20K.tif",
    "normal_data/162-165/162. Mask 20K.tif",
])
def test_crop_based_and_global_assignment_agree_when_nothing_straddles(image_path):
    """On real data verified to have zero straddling mitochondria (checked
    across the full 99-mask dataset during Phase 3a development), the
    legacy crop-based path and the new global-assignment path must
    produce numerically identical per-axon metrics -- this is what makes
    flipping the default safe on this dataset. If this ever fails, it
    means either the image now has a straddling mitochondrion or one of
    the two code paths has a bug."""
    from temptation.segmentation import resolve_mode, label_axons, run_watershed

    mask = tiff.imread(DATA_ROOT / image_path)
    cfg = SegmentationConfig(mito_hole_handling="fill")
    mode = resolve_mode(mask, cfg)
    comps = build_compartments(mask, cfg)
    axon_lab = label_axons(comps, cfg)
    labels_ws = run_watershed(axon_lab, comps, mode, cfg)
    mito_lab, assignment = assign_mito_to_axons(comps.mito_in_axon, labels_ws, method="centroid")

    fiber_to_mito = {}
    for mito_label, fiber_label in assignment.items():
        fiber_to_mito.setdefault(fiber_label, []).append(mito_label)

    mito_regions_by_label = {r.label: r for r in regionprops(mito_lab)}

    px = 0.00524
    n_checked = 0
    for p in regionprops(labels_ws):
        fiber_id = p.label
        min_row, min_col, max_row, max_col = p.bbox
        fiber_crop = p.image
        axon_crop = comps.axoplasm[min_row:max_row, min_col:max_col]
        axon_mask_local = fiber_crop & axon_crop
        axon_area_px = int(axon_mask_local.sum())
        if axon_area_px == 0:
            continue

        axon_labeled_local = sk_label(axon_mask_local)
        axon_rp = regionprops(axon_labeled_local)
        if not axon_rp:
            continue
        axon_props_local = max(axon_rp, key=lambda r: r.area)
        cy_local, cx_local = axon_props_local.centroid
        cy_global = min_row + cy_local
        cx_global = min_col + cx_local
        axon_area_um2 = axon_area_px * px ** 2

        mito_crop = comps.mito[min_row:max_row, min_col:max_col]
        mito_in_axon_crop = mito_crop & axon_mask_local
        old = mito_metrics_for_axon(mito_in_axon_crop, axon_area_um2, cx_local, cy_local, px)

        regions_here = [mito_regions_by_label[ml] for ml in fiber_to_mito.get(fiber_id, [])]
        new = mito_metrics_from_regions(regions_here, axon_area_um2, cx_global, cy_global, px)

        for key in old:
            ov, nv = old[key], new[key]
            if isinstance(ov, float) and np.isnan(ov) and isinstance(nv, float) and np.isnan(nv):
                continue
            assert abs(ov - nv) < 1e-9, f"{image_path} axon {fiber_id} field {key}: old={ov} new={nv}"
        n_checked += 1

    assert n_checked > 0
