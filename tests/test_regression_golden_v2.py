"""Regression anchor for the *corrected* (fill-mode, default since Phase
2b) output. Companion to test_regression_golden.py, which pins the
pre-fix numbers under --mito-hole-handling legacy. If this test fails, a
later phase has changed the fixed F1 behavior without updating this
golden set deliberately -- do not "fix" this test by just re-copying
current output over tests/golden_v2/ without checking why the numbers
moved."""

from pathlib import Path

import pandas as pd
import pytest

from temptation.cli import main as cli_main

GOLDEN_V2 = Path(__file__).parent / "golden_v2"
DATA_ROOT = Path(__file__).parent.parent.parent

GROUP_ARGS = {
    "normal": ["--folder", str(DATA_ROOT / "normal_data" / "162-165")],
    "pathological": ["--folder", str(DATA_ROOT / "pathological_data")],
}


# image_path/mask_path record whatever path was passed on the command
# line (relative or absolute depending on how the CLI was invoked) --
# that's invocation-dependent metadata, not a computed value, so it is
# deliberately excluded from the byte-equality check below.
_PATH_COLUMNS = ("image_path", "mask_path")


@pytest.mark.parametrize("group", ["normal", "pathological"])
def test_axons_csv_matches_golden_v2(group, tmp_path, monkeypatch):
    out_dir = tmp_path / group
    out_dir.mkdir()
    argv = ["prog", *GROUP_ARGS[group], "--pixel-um", "0.00524", "--output-dir", str(out_dir)]
    monkeypatch.setattr("sys.argv", argv)
    cli_main()

    new = pd.read_csv(out_dir / "axons.csv")
    golden = pd.read_csv(GOLDEN_V2 / group / "axons.csv")

    assert list(new.columns) == list(golden.columns)
    value_cols = [c for c in new.columns if c not in _PATH_COLUMNS]
    pd.testing.assert_frame_equal(new[value_cols], golden[value_cols], check_exact=False, rtol=1e-12)
    for c in _PATH_COLUMNS:
        assert new[c].notna().all()


@pytest.mark.parametrize("group", ["normal", "pathological"])
def test_image_summary_csv_matches_golden_v2(group, tmp_path, monkeypatch):
    out_dir = tmp_path / group
    out_dir.mkdir()
    argv = ["prog", *GROUP_ARGS[group], "--pixel-um", "0.00524", "--output-dir", str(out_dir)]
    monkeypatch.setattr("sys.argv", argv)
    cli_main()

    new = pd.read_csv(out_dir / "image_summary.csv")
    golden = pd.read_csv(GOLDEN_V2 / group / "image_summary.csv")

    assert list(new.columns) == list(golden.columns)
    value_cols = [c for c in new.columns if c not in _PATH_COLUMNS]
    pd.testing.assert_frame_equal(new[value_cols], golden[value_cols], check_exact=False, rtol=1e-12)
    for c in _PATH_COLUMNS:
        assert new[c].notna().all()


@pytest.mark.parametrize("group", ["normal", "pathological"])
def test_fill_mode_axon_counts_match_legacy_mode(group, tmp_path, monkeypatch):
    """On these two specific folders (the golden-baseline regression
    anchor), axon counts happen to be identical between modes. This is
    NOT a general guarantee -- see test_fill_mode_can_merge_oversegmented_
    axons below and docs/phase2b_impact.md for 7/100 real-dataset images
    where fill mode legitimately produces one fewer axon (a large
    mitochondrion bisecting an axon's footprint under the legacy bug was
    causing watershed to over-segment it into two pieces). Pinned here so
    a future unrelated change to masks.py that perturbs axon_lab seeding
    on *these* two folders specifically is caught."""
    out_dir = tmp_path / group
    out_dir.mkdir()
    argv = ["prog", *GROUP_ARGS[group], "--pixel-um", "0.00524", "--output-dir", str(out_dir)]
    monkeypatch.setattr("sys.argv", argv)
    cli_main()

    fill_axons = pd.read_csv(out_dir / "axons.csv")
    legacy_axons = pd.read_csv(Path(__file__).parent / "golden" / group / "axons.csv")
    assert len(fill_axons) == len(legacy_axons)
    assert sorted(fill_axons.image_id.astype(str)) == sorted(legacy_axons.image_id.astype(str))


def test_fill_mode_can_merge_oversegmented_axons(tmp_path, monkeypatch):
    """Documents and pins the real (dataset-wide, not golden-anchor-local)
    exception to the above: image 172 in normal_data/171-175 has 16 axons
    under legacy mode and 15 under fill mode. Verified mechanism (see
    docs/phase2b_impact.md): under legacy, a large mitochondrion bisects
    one axon's axon_only footprint, so axoplasm is disconnected there and
    watershed produces a real fiber plus a spurious ~900px fragment; under
    fill, that same mitochondrion (99.2% of the connecting gap) is
    correctly recognized as part of the axon and the two pieces merge
    into one properly-sized fiber. This is the expected, beneficial
    behavior of the fix -- if this test starts failing because the count
    difference disappears, someone should check whether the merge
    mechanism still fires as documented, not just update the numbers."""
    image_dir = DATA_ROOT / "normal_data" / "171-175"

    out_fill = tmp_path / "fill"
    out_fill.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(image_dir), "--pixel-um", "0.00524",
        "--output-dir", str(out_fill),
    ])
    cli_main()
    fill_axons = pd.read_csv(out_fill / "axons.csv")

    out_legacy = tmp_path / "legacy"
    out_legacy.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(image_dir), "--pixel-um", "0.00524",
        "--output-dir", str(out_legacy), "--mito-hole-handling", "legacy",
    ])
    cli_main()
    legacy_axons = pd.read_csv(out_legacy / "axons.csv")

    n_fill_172 = len(fill_axons[fill_axons.image_id == 172])
    n_legacy_172 = len(legacy_axons[legacy_axons.image_id == 172])
    assert n_legacy_172 == 16
    assert n_fill_172 == 15
    assert n_legacy_172 - n_fill_172 == 1


def test_schema_version_reflects_active_mito_hole_handling(tmp_path, monkeypatch):
    from temptation.schema import SCHEMA_VERSION, SCHEMA_VERSION_FILL

    out_fill = tmp_path / "fill"
    out_fill.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_fill),
    ])
    cli_main()
    df_fill = pd.read_csv(out_fill / "axons.csv")
    assert (df_fill["schema_version"].astype(str) == SCHEMA_VERSION_FILL).all()

    out_legacy = tmp_path / "legacy"
    out_legacy.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_legacy),
        "--mito-hole-handling", "legacy",
    ])
    cli_main()
    df_legacy = pd.read_csv(out_legacy / "axons.csv")
    assert (df_legacy["schema_version"].astype(str) == SCHEMA_VERSION).all()
