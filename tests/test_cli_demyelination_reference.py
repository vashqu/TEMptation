"""CLI integration test for --demyelination-reference (Phase 4)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from temptation.cli import main as cli_main

DATA_ROOT = Path(__file__).parent.parent.parent


def test_demyelination_index_nan_by_default(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
    ])
    cli_main()
    df = pd.read_csv(out_dir / "image_summary.csv")
    assert df["image_demyelination_index"].isna().all()


def test_demyelination_index_computed_when_reference_supplied(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
        "--demyelination-reference", "0.30",
    ])
    cli_main()
    df = pd.read_csv(out_dir / "image_summary.csv")
    assert df["image_demyelination_index"].notna().all()
    assert (df["image_demyelination_index"] >= 0).all()
    assert (df["image_demyelination_index"] <= 1).all()

    expected = (1.0 - df["myelin_area_fraction_of_fov"] / 0.30).clip(0.0, 1.0)
    pd.testing.assert_series_equal(
        df["image_demyelination_index"], expected, check_names=False, rtol=1e-9,
    )
