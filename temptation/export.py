"""All CSV/JSON writes go through this module -- called identically by
CLI and GUI so file naming never drifts between the two front-ends."""

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def _warn_if_exists(path: Path) -> None:
    """Overwrite guard (Phase 6): not a hard block (this tool has no
    interactive prompt to block on, and batch/automation use must not
    silently fail), but every overwrite is announced so it's visible in
    the run's console output."""
    if path.exists():
        print(f"  [WARNING] Overwriting existing file: {path}")


def write_axon_csv(df_axons: pd.DataFrame, out_dir: Path) -> Path:
    path = Path(out_dir) / "axons.csv"
    _warn_if_exists(path)
    df_axons.to_csv(path, index=False)
    return path


def write_image_csv(df_images: pd.DataFrame, out_dir: Path) -> Path:
    path = Path(out_dir) / "image_summary.csv"
    _warn_if_exists(path)
    df_images.to_csv(path, index=False)
    return path


def write_mito_csv(df_mito: pd.DataFrame, out_dir: Path) -> Path:
    """One row per real mitochondrion (Phase 3b, --write-mito-csv)."""
    path = Path(out_dir) / "mitochondria_metrics.csv"
    _warn_if_exists(path)
    df_mito.to_csv(path, index=False)
    return path


def write_qc_report_csv(df_qc: pd.DataFrame, out_dir: Path) -> Path:
    """Slim per-axon QC view (Phase 5, --write-qc-report)."""
    path = Path(out_dir) / "qc_report.csv"
    _warn_if_exists(path)
    df_qc.to_csv(path, index=False)
    return path


def write_group_csv(df_groups: pd.DataFrame, out_dir: Path) -> Path:
    """One row per group (Phase 6, --write-group-csv). Built from
    summaries.summarize_groups()."""
    path = Path(out_dir) / "group_metrics.csv"
    _warn_if_exists(path)
    df_groups.to_csv(path, index=False)
    return path


def _package_versions() -> dict:
    versions = {}
    for display_name, import_name in (
        ("numpy", "numpy"), ("scipy", "scipy"), ("scikit-image", "skimage"),
        ("pandas", "pandas"), ("tifffile", "tifffile"), ("matplotlib", "matplotlib"),
    ):
        try:
            mod = __import__(import_name)
            versions[display_name] = getattr(mod, "__version__", "unknown")
        except ImportError:
            versions[display_name] = "not installed"
    return versions


def write_run_manifest(
    out_dir: Path,
    config: dict,
    inputs: list,
    skipped_files: list,
    outputs: dict,
    row_counts: dict,
    exclusion_counts: dict,
    processing_errors: list = None,
) -> Path:
    """run_manifest.json (CLAUDE.md Sec 16.8, Phase 6): everything needed
    to reproduce or audit a run without re-reading console output --
    config (every CLI argument), package versions, a UTC timestamp, every
    input pair actually processed, every file skipped during discovery,
    which output files were written and how many rows each has, a
    per-rule breakdown of why axons were excluded (qc.exclusion_reason_counts),
    and any per-image processing errors (pipeline.DatasetResult.errors --
    a failed image doesn't abort the batch, so this is the only durable
    record of which images silently contributed nothing).

    All values are passed in already-computed rather than recomputed here
    -- this function only serializes; it never re-derives anything, so it
    can't silently disagree with what was actually written to the CSVs.
    """
    from . import __version__ as temptation_version

    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "temptation_version": temptation_version,
        "package_versions": _package_versions(),
        "config": config,
        "inputs": inputs,
        "skipped_files": skipped_files,
        "outputs": {k: str(v) for k, v in outputs.items()},
        "row_counts": row_counts,
        "exclusion_counts": exclusion_counts,
        "processing_errors": processing_errors or [],
    }
    path = Path(out_dir) / "run_manifest.json"
    _warn_if_exists(path)
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    return path
