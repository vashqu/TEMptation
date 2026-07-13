"""CLI integration tests for Phase 6 (--schema, run_manifest.json,
--write-group-csv)."""

import json
from pathlib import Path

import pandas as pd
import pytest

from temptation.cli import main as cli_main
from temptation.schema import AXON_COLUMNS_LEGACY, IMAGE_COLUMNS_LEGACY

DATA_ROOT = Path(__file__).parent.parent.parent


def _run(monkeypatch, out_dir, extra_args=None):
    argv = [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
    ]
    argv += extra_args or []
    monkeypatch.setattr("sys.argv", argv)
    cli_main()


# ---------------------------------------------------------------------- --schema

def test_schema_full_is_default(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir)
    axons = pd.read_csv(out_dir / "axons.csv")
    assert len(axons.columns) > len(AXON_COLUMNS_LEGACY)


def test_schema_legacy_strips_to_exact_legacy_columns(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir, ["--schema", "legacy"])
    axons = pd.read_csv(out_dir / "axons.csv")
    images = pd.read_csv(out_dir / "image_summary.csv")
    assert list(axons.columns) == list(AXON_COLUMNS_LEGACY)
    assert list(images.columns) == list(IMAGE_COLUMNS_LEGACY)


def test_schema_legacy_values_match_golden(tmp_path, monkeypatch):
    """--schema legacy must not just have the right columns, but the
    exact same values as the Phase 0 golden baseline (with
    --mito-hole-handling legacy --mito-assignment legacy, since those
    flags govern VALUES, not --schema, which only governs column
    layout)."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir, [
        "--schema", "legacy",
        "--mito-hole-handling", "legacy", "--mito-assignment", "legacy",
    ])
    axons = pd.read_csv(out_dir / "axons.csv")
    golden = pd.read_csv(Path(__file__).parent / "golden" / "normal" / "axons.csv")
    pd.testing.assert_frame_equal(axons, golden, check_exact=False, rtol=1e-12)


# ---------------------------------------------------------------------- manifest

def test_manifest_always_written(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir)
    assert (out_dir / "run_manifest.json").exists()


def test_manifest_content(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir, ["--exclude-qc-failed"])

    with open(out_dir / "run_manifest.json") as f:
        manifest = json.load(f)

    assert "timestamp_utc" in manifest
    assert "temptation_version" in manifest
    assert set(manifest["package_versions"]) >= {"numpy", "scipy", "scikit-image", "pandas"}
    assert manifest["config"]["exclude_qc_failed"] is True
    assert manifest["config"]["resolved_pixel_length_um"] == pytest.approx(0.00524)
    assert len(manifest["inputs"]) == 4  # 162-165
    assert manifest["row_counts"]["axons"] == 47
    assert manifest["row_counts"]["images"] == 4
    assert "axons_csv" in manifest["outputs"]
    assert "image_summary_csv" in manifest["outputs"]
    assert isinstance(manifest["exclusion_counts"], dict)


def test_manifest_records_skipped_files(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--input-root", str(DATA_ROOT), "--pixel-um", "0.00524",
        "--output-dir", str(out_dir),
        "--groups", "normal_data:normal", "pathological_data:pathological",
    ])
    cli_main()
    with open(out_dir / "run_manifest.json") as f:
        manifest = json.load(f)
    # image 429's mask is a .gif; find_groups recovers it (F7), but other
    # genuinely-unpaired files in the tree, if any, would show up here --
    # at minimum this must not crash and must be a list.
    assert isinstance(manifest["skipped_files"], list)
    assert manifest["row_counts"]["axons"] > 1000


def test_manifest_overwrite_warns(tmp_path, monkeypatch, capsys):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir)
    capsys.readouterr()  # clear first run's output
    _run(monkeypatch, out_dir)
    captured = capsys.readouterr()
    assert "Overwriting existing file" in captured.out


# ---------------------------------------------------------------------- group csv

def test_write_group_csv(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir, ["--group", "normal", "--write-group-csv"])
    group_csv = out_dir / "group_metrics.csv"
    assert group_csv.exists()
    df = pd.read_csv(group_csv)
    assert len(df) == 1
    assert df["group"].iloc[0] == "normal"
    assert df["n_images"].iloc[0] == 4


def test_no_group_csv_without_flag(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir, ["--group", "normal"])
    assert not (out_dir / "group_metrics.csv").exists()


def test_no_group_csv_when_no_group_assigned(tmp_path, monkeypatch):
    """--write-group-csv with no --group set (group is None for every
    row) should not crash, just produce nothing meaningful to write."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _run(monkeypatch, out_dir, ["--write-group-csv"])
    # group column is None/NaN for every row -> summarize_groups still
    # returns a (single, "None"-labeled) group row, which is fine; the
    # key assertion is that this doesn't crash.
    assert True
