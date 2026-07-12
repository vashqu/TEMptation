"""Tests for --write-mito-csv (mitochondria_metrics.csv, Phase 3b)."""

from pathlib import Path

import pandas as pd
import pytest

from temptation.cli import main as cli_main
from temptation.schema import MITO_COLUMNS

DATA_ROOT = Path(__file__).parent.parent.parent


def test_write_mito_csv_default_assignment(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir), "--write-mito-csv",
    ])
    cli_main()

    mito_csv = out_dir / "mitochondria_metrics.csv"
    assert mito_csv.exists()
    df = pd.read_csv(mito_csv)
    assert list(df.columns)[:len(MITO_COLUMNS)] == list(MITO_COLUMNS) or set(MITO_COLUMNS).issubset(df.columns)
    assert len(df) > 0

    # one row per mitochondrion, not per fragment: total rows should equal
    # the sum of mito_count across axons.csv for images in this folder
    axons = pd.read_csv(out_dir / "axons.csv")
    assert len(df) == int(axons["mito_count"].sum())

    # every mito row's parent_axon_id must correspond to a real axon_id
    # within the same image
    merged = df.merge(
        axons[["image_id", "axon_id"]].drop_duplicates(),
        left_on=["image_id", "parent_axon_id"], right_on=["image_id", "axon_id"],
        how="left",
    )
    assert merged["axon_id"].notna().all()


def test_write_mito_csv_empty_under_legacy_assignment(tmp_path, monkeypatch, capsys):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
        "--write-mito-csv", "--mito-assignment", "legacy",
    ])
    cli_main()
    mito_csv = out_dir / "mitochondria_metrics.csv"
    assert not mito_csv.exists()
    captured = capsys.readouterr()
    assert "legacy" in captured.out.lower()


def test_no_mito_csv_written_without_flag(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
    ])
    cli_main()
    assert not (out_dir / "mitochondria_metrics.csv").exists()


def test_mito_csv_area_sum_matches_axon_level_total(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir), "--write-mito-csv",
    ])
    cli_main()
    mito = pd.read_csv(out_dir / "mitochondria_metrics.csv")
    axons = pd.read_csv(out_dir / "axons.csv")

    per_axon_mito_sum = mito.groupby(["image_id", "parent_axon_id"])["mito_area_um2"].sum()
    for (image_id, axon_id), total in per_axon_mito_sum.items():
        expected = axons[(axons.image_id == image_id) & (axons.axon_id == axon_id)]["mito_total_area_um2"].iloc[0]
        assert total == pytest.approx(expected, rel=1e-9)
