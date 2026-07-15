# Changelog

All notable changes to TEMptation, in [Keep a Changelog](https://keepachangelog.com/)
style. Dates are omitted (this project doesn't do dated releases); entries are grouped
by development phase and appear oldest-first.

**If you're deciding whether to trust a number you already exported**: check that
row's `schema_version` column. `"2.0"` means the Phase 2b mitochondrial-area fix was
active; `"1.0"` means it wasn't. See `docs/known_issues.md`.

## Pre-refactor baseline

The original tool: a single-file CLI (`measure_nerve.py`) and a single-file GUI
(`measure_nerve_gui.py`), sharing one `measure_image()` function. Frozen as
`tests/golden/` before any further changes, and re-verified in `AUDIT.md` — eight
defects (F1–F8) documented there, most fixed in later phases below.

## Phase 1 — Backend package extraction

Split the monolithic script into `temptation/` (`config`, `mathutils`, `dataio`,
`discovery`, `masks`, `segmentation`, `metrics_axon`, `metrics_mito`, `metrics_spatial`,
`schema`, `export`, `plotting`, `pipeline`, `cli`, `compat`). Zero output change —
every function moved verbatim; `axons.csv`/`image_summary.csv` byte-identical to the
frozen golden baseline.

## Phase 2 — Metadata, recursive discovery, unbiased shape metrics

- Added identity/metadata columns: `group`, `image_path`, `mask_path`, `pixel_size_um`,
  `schema_version`.
- `--input-root`/`--groups`/`--group`/`--recursive` for multi-group batch runs.
- Recursive folder scanning with an extended glob (`.tif`/`.tiff`/`.png`/`.gif`) that
  reports skipped files instead of silently dropping them (fixes F7).
- `axon_circularity_crofton`, `axon_shape_irregularity`, `axon_shape_irregularity_crofton`
  — an unbiased alternative to the legacy `circularity` column, added alongside it
  (fixes F5's documentation gap; the legacy column itself is unchanged).

## Phase 2b — Mitochondrial-rim bug fix (gated, values change)

**The only phase that changes previously-exported numbers.** Fixed F1 (mitochondrial
area was ~2% of true value) and F2 (`mito_outside_frac` structurally always `1.0`) via
`--mito-hole-handling {fill,legacy}`, `fill` now the default. `legacy` reproduces the
original bug exactly, for regression testing only. `schema_version` bumps to `"2.0"`
when `fill` is active. See `docs/known_issues.md` F1/F2 and `docs/phase2b_impact.md`
for the full before/after numbers.

## Phase 3a — Mitochondrion double-counting fix

Fixed F8: mitochondria straddling a watershed boundary were fragmented and
double-counted by the original per-fiber-crop approach. `--mito-assignment
{centroid,overlap,legacy}`, `centroid` now the default — mitochondria are labeled once
globally, then each whole mitochondrion is assigned to exactly one axon.
`n_mito_assigned`/`n_mito_unassigned` diagnostics added to `image_summary.csv`.

## Phase 3b — Mitochondrial burden, shape, and spatial metrics

Added `mito_total_area_um2`, `mito_occupancy_ratio`, `mito_mean_area_um2`,
`mito_median_area_um2`, `mito_area_iqr`, `mito_fragmentation_index`,
`normalized_mito_load`, `mito_per_myelin`, `mito_mean_aspect_ratio`,
`mito_std_aspect_ratio`, `mito_mean_solidity`, `mito_std_solidity`,
`mito_mean_eccentricity`, `mito_peripheralization_index`, `mito_mean_nn_distance_um`,
`mito_clustering_index`. New `mitochondria_metrics.csv` output (`--write-mito-csv`,
one row per real mitochondrion — requires `--mito-assignment centroid` or `overlap`).
Full column reference: `docs/metrics.md`.

## Phase 4 — Image-level distribution statistics

Added `image_{mean,std,cv,median,min,max,iqr}_{g_ratio,axon_area_um2,
mito_density_per_um2,mito_occupancy_ratio}` (28 columns), `image_mvf`, and
`image_demyelination_index` (`--demyelination-reference`, no default — a biological
calibration choice this tool cannot supply on its own; `NaN` until you provide one).
See `docs/known_issues.md` F4 for why a naive demyelination index fails on this
tool's own pathological test data.

## Phase 5 — QC flags and exclusion marking

Added every `qc_*` flag (always computed), `excluded_from_analysis`/`exclusion_reason`
(only populated when `--exclude-qc-failed` is set; rows are marked, never dropped),
and `qc_report.csv` (`--write-qc-report`). Full reference: `docs/qc.md`, including the
three flags (`qc_no_myelin`, `qc_no_mitochondria`, `qc_low_mito_count_clustering`)
that can never trigger exclusion by design.

## Phase 6 — Schema stability and export redesign

Fixed the other half of F8: every output is now reindexed to a fixed column set
(`--schema {legacy,full}`) regardless of how many objects were detected — a 0-axon
image emits the full header, not a truncated one. Added `run_manifest.json` (always
written: config, package versions, timestamp, inputs, skipped files, outputs, row
counts, exclusion counts) and `group_metrics.csv` (`--write-group-csv`).

## Phase 7 — GUI modernization

Replaced the single-window GUI with a step-rail application (`gui/`), built as a thin
interface layer over `temptation.pipeline` — enforced by an `ast`-based test
(`tests/test_gui_has_no_metrics.py`) that fails if any metric-computing call appears
in GUI code. Landed across five sub-phases: shared batch orchestration
(`pipeline.analyze_dataset`), the dataset/calibration/run shell, metric-selection and
QC panels, the results dashboard and export panel, and the visual review (overlay +
per-axon inspector) panel. The pre-refactor GUI is retained as
`measure_nerve_gui_legacy.py`.

## Post-launch GUI revision — navigation, speed, clarity

Real usage surfaced issues the ast-guard and unit tests couldn't catch:

- **Navigation reorganized**: 7 rail steps collapsed to 3 (Setup/Review/Export).
  Setup bundles Data/Calibration/Segmentation/QC as scrollable sub-tabs with Run
  pinned below; output directory moved to Export (Run never wrote a file — only
  Export does).
- **A real responsiveness bug fixed**, not just perceived slowness: editing a QC
  threshold used to null the last run's results on every keystroke, blanking the
  dashboard and export checklist mid-typing. Now edits set a `result_stale` flag
  instead — results stay visible with a "parameters changed since this run" note.
  Input handlers debounced; dashboard/export/QC-preview rebuilds gated on whether
  their actual dependency changed, not on every unrelated edit; mask-availability
  scanning moved to a background thread; the Inspect panel's per-image QC badges
  now use one vectorized pass instead of one dataframe filter per image.
- **Inspect panel layout bug fixed**: the embedded overlay figure was sized for its
  own full popup window, crowding the image list and axon inspector down to
  near-zero width once an image loaded. Fixed with a smaller embedded figure size,
  protected column widths, and explicit Previous/Next buttons.
- **`mode` vs `group` clarified**: a labeled "Auto-detect / Control / Pathology"
  selector replaces a raw parameter-accordion entry, and the `mode` column is now
  optional in GUI exports (off by default). See `docs/known_issues.md` F3.
- Dashboard summary-card layout gap fixed (a row-weighting conflict between two
  files, not the dashboard logic itself).

## Not yet implemented

- **Morphometric pathology score** — planned, not built. See `docs/pathology_score.md`.
- **`--config`/config-file loading** — `temptation.config.load_config()` exists but
  isn't wired to the CLI or GUI. See `config.example.yaml`'s header.
- **Pixel-space-only export mode** (CLAUDE.md §16.3: `_px`-suffixed columns, no
  `_um2` columns, for uncalibrated data) — the GUI's "pixel-space only" checkbox is a
  documented placeholder today; it runs with an identity scale rather than actually
  restructuring output columns.
- **Biological validation on real data** (Phase 9) — not yet run. See
  `docs/pathology_score.md`'s note on what "done" looks like for that phase, and the
  same discipline applies to any future validation work: descriptive checks only, no
  hypothesis tests, no classifier accuracy figures.
