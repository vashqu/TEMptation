"""Acceptance test: the refactored/extended package must still reproduce
the pre-refactor golden CSVs on the *legacy* columns exactly (Phase 2 adds
new columns -- axon shape metrics, group/path/pixel_size/schema_version
identity metadata -- so a plain equality check on the full column set no
longer applies; see IMPLEMENTATION_BLUEPRINT.md Sec 9, Phase 2 "Checks").
If the legacy-subset comparison below fails, a value that used to be
exported has silently changed -- that is not allowed outside the gated
Phase 2b."""

from pathlib import Path

import pandas as pd
import pytest

from temptation.cli import main as cli_main

GOLDEN = Path(__file__).parent / "golden"
DATA_ROOT = Path(__file__).parent.parent.parent

GROUP_ARGS = {
    "normal": ["--folder", str(DATA_ROOT / "normal_data" / "162-165")],
    "pathological": ["--folder", str(DATA_ROOT / "pathological_data")],
}


@pytest.mark.parametrize("group", ["normal", "pathological"])
def test_axons_csv_matches_golden(group, tmp_path, monkeypatch):
    out_dir = tmp_path / group
    out_dir.mkdir()
    argv = ["prog", *GROUP_ARGS[group], "--pixel-um", "0.00524", "--output-dir", str(out_dir)]
    monkeypatch.setattr("sys.argv", argv)
    cli_main()

    new = pd.read_csv(out_dir / "axons.csv")
    golden = pd.read_csv(GOLDEN / group / "axons.csv")

    assert set(golden.columns).issubset(set(new.columns)), (
        f"missing legacy columns: {set(golden.columns) - set(new.columns)}"
    )
    pd.testing.assert_frame_equal(new[golden.columns], golden, check_exact=False, rtol=1e-12)


@pytest.mark.parametrize("group", ["normal", "pathological"])
def test_image_summary_csv_matches_golden(group, tmp_path, monkeypatch):
    out_dir = tmp_path / group
    out_dir.mkdir()
    argv = ["prog", *GROUP_ARGS[group], "--pixel-um", "0.00524", "--output-dir", str(out_dir)]
    monkeypatch.setattr("sys.argv", argv)
    cli_main()

    new = pd.read_csv(out_dir / "image_summary.csv")
    golden = pd.read_csv(GOLDEN / group / "image_summary.csv")

    assert set(golden.columns).issubset(set(new.columns)), (
        f"missing legacy columns: {set(golden.columns) - set(new.columns)}"
    )
    pd.testing.assert_frame_equal(new[golden.columns], golden, check_exact=False, rtol=1e-12)


@pytest.mark.parametrize("group,n_axon_rows,n_image_rows", [
    ("normal", 47, 4),
    ("pathological", 25, 6),
])
def test_row_counts_match_baseline(group, n_axon_rows, n_image_rows, tmp_path, monkeypatch):
    out_dir = tmp_path / group
    out_dir.mkdir()
    argv = ["prog", *GROUP_ARGS[group], "--pixel-um", "0.00524", "--output-dir", str(out_dir)]
    monkeypatch.setattr("sys.argv", argv)
    cli_main()

    axons = pd.read_csv(out_dir / "axons.csv")
    images = pd.read_csv(out_dir / "image_summary.csv")
    assert len(axons) == n_axon_rows
    assert len(images) == n_image_rows
