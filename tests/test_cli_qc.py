"""CLI integration tests for Phase 5 QC/exclusion flags and qc_report.csv."""

from pathlib import Path

import pandas as pd
import pytest

from temptation.cli import main as cli_main

DATA_ROOT = Path(__file__).parent.parent.parent


def test_default_run_excludes_nothing(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
    ])
    cli_main()
    df = pd.read_csv(out_dir / "axons.csv")
    assert (df["excluded_from_analysis"] == False).all()  # noqa: E712
    assert (df["exclusion_reason"].fillna("") == "").all()


def test_row_count_invariant_to_exclude_flag(tmp_path, monkeypatch):
    counts = {}
    for exclude in (False, True):
        out_dir = tmp_path / f"out_{exclude}"
        out_dir.mkdir()
        argv = [
            "prog", "--folder", str(DATA_ROOT / "pathological_data"),
            "--pixel-um", "0.00524", "--output-dir", str(out_dir),
        ]
        if exclude:
            argv.append("--exclude-qc-failed")
        monkeypatch.setattr("sys.argv", argv)
        cli_main()
        counts[exclude] = len(pd.read_csv(out_dir / "axons.csv"))
    assert counts[False] == counts[True]


def test_pathological_group_retains_valid_axons_with_exclusion(tmp_path, monkeypatch):
    """The Sec 6.5 safety property, at the CLI level: every pathological
    image must retain at least one valid axon, and qc_no_myelin/
    invalid_area_relation-from-zero-myelin must never appear in
    exclusion_reason."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "pathological_data"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
        "--exclude-qc-failed",
    ])
    cli_main()
    images = pd.read_csv(out_dir / "image_summary.csv")
    assert (images["image_n_axons_valid"] > 0).all()

    axons = pd.read_csv(out_dir / "axons.csv")
    reasons = axons.loc[axons["excluded_from_analysis"], "exclusion_reason"].dropna()
    assert not reasons.str.contains("myelin").any()


def test_write_qc_report(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
        "--write-qc-report",
    ])
    cli_main()
    qc_csv = out_dir / "qc_report.csv"
    assert qc_csv.exists()
    report = pd.read_csv(qc_csv)
    axons = pd.read_csv(out_dir / "axons.csv")
    assert len(report) == len(axons)
    assert "qc_g_ratio_low" in report.columns
    assert "mito_count" not in report.columns  # slim view, not the full metric set


def test_no_qc_report_without_flag(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
        "--pixel-um", "0.00524", "--output-dir", str(out_dir),
    ])
    cli_main()
    assert not (out_dir / "qc_report.csv").exists()


def test_custom_qc_thresholds(tmp_path, monkeypatch):
    """A very permissive circularity threshold should flag more axons
    than the strict default."""
    def run(threshold):
        out_dir = tmp_path / f"out_{threshold}"
        out_dir.mkdir()
        monkeypatch.setattr("sys.argv", [
            "prog", "--folder", str(DATA_ROOT / "normal_data" / "162-165"),
            "--pixel-um", "0.00524", "--output-dir", str(out_dir),
            "--axon-circularity-min", str(threshold),
        ])
        cli_main()
        df = pd.read_csv(out_dir / "axons.csv")
        return df["qc_low_circularity"].sum()

    n_strict = run(0.0)   # nothing can be below 0.0
    n_loose = run(1.0)    # everything is below 1.0
    assert n_strict == 0
    assert n_loose > n_strict
