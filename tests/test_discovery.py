"""Tests for the Phase 2 discovery additions (F7 fixes:
IMPLEMENTATION_BLUEPRINT.md Sec 0). `find_pairs` (Phase 1, unchanged) is
covered by test_regression_golden.py; this file covers find_groups,
infer_group, and parse_group_map."""

from pathlib import Path

import pytest

from temptation.discovery import find_groups, find_pairs, infer_group, parse_group_map

DATA_ROOT = Path(__file__).parent.parent.parent


def test_parse_group_map():
    gm = parse_group_map(["normal_data:normal", "pathological_data:pathological"])
    assert gm == {"normal_data": "normal", "pathological_data": "pathological"}


def test_parse_group_map_rejects_missing_colon():
    with pytest.raises(ValueError):
        parse_group_map(["normal_data"])


def test_find_groups_recovers_recursive_pairs_and_extended_glob():
    group_map = {"normal_data": "normal", "pathological_data": "pathological"}
    pairs, skipped = find_groups(DATA_ROOT, group_map)

    n_normal = sum(1 for p in pairs if p["group"] == "normal")
    n_path = sum(1 for p in pairs if p["group"] == "pathological")

    assert n_normal == 94
    assert n_path == 6

    ids = {(p["id"], p["group"]) for p in pairs}
    assert len(ids) == len(pairs), "no duplicate (id, group) pairs"

    # F7: image 429's mask is a .gif and used to be silently dropped by the
    # non-recursive, *.tif/*.tiff-only find_pairs. The extended glob must
    # recover it as a real pair, not report it as skipped.
    recovered_429 = [p for p in pairs if p["group"] == "normal" and p["id"] == "429"]
    assert len(recovered_429) == 1
    assert recovered_429[0]["mask_path"].suffix == ".gif"
    gif_skips = [s for s in skipped if s.path.suffix == ".gif"]
    assert gif_skips == []


def test_find_groups_reports_missing_folder_as_skipped_not_an_exception():
    pairs, skipped = find_groups(DATA_ROOT, {"does_not_exist_data": "nowhere"})
    assert pairs == []
    assert len(skipped) == 1
    assert "not found" in skipped[0].reason


def test_infer_group():
    group_map = {"normal_data": "normal", "pathological_data": "pathological"}
    assert infer_group(DATA_ROOT / "normal_data" / "162-165" / "163. Mask 20K.tif", group_map) == "normal"
    assert infer_group(DATA_ROOT / "pathological_data" / "mask_path_1.tif", group_map) == "pathological"
    assert infer_group(Path("/some/unrelated/path.tif"), group_map) is None


def test_find_pairs_phase1_behavior_is_unaffected_by_phase2_additions():
    """find_pairs must still drop 429 (non-recursive, *.tif/*.tiff-only) --
    this is the exact frozen Phase 1 behavior that --folder (no
    --recursive) still relies on."""
    pairs = find_pairs(DATA_ROOT / "normal_data" / "422-430")
    assert len(pairs) == 8
    assert not any(p["id"] == "429" for p in pairs)
