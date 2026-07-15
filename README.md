# TEMptation — Nerve Morphometry

Segment and measure nerve fibers from **TEM images + segmentation masks**.
Extracts axon geometry, myelin/g-ratio metrics, mitochondrial burden/shape/spatial
metrics, image-level distribution statistics, and QC flags — for **normal**,
**pathological**, and (via a labeled `group`) any other experimental condition.

Two interfaces, one engine: a **CLI** (`measure_nerve.py`) and a **GUI**
(`measure_nerve_gui.py`), both calling the same `temptation/` package. Neither
reimplements the other's logic — see `docs/known_issues.md` and
`tests/test_gui_has_no_metrics.py` if you're curious how that's enforced.

> **Not a diagnostic tool.** Every metric here is a morphometric measurement derived
> from a segmentation mask. None of it is validated for, or intended for, clinical
> or research diagnosis. See `docs/pathology_score.md` for more on this.

---

## Documentation map

| Doc | What's in it |
|---|---|
| This file | Install, quick start, CLI/GUI overview, troubleshooting |
| `how_to_use.txt` | The same content in longer, plain-text, step-by-step form |
| `docs/metrics.md` | Every output column: formula, interpretation, NaN conditions |
| `docs/qc.md` | Every QC flag, thresholds, what can and can't trigger exclusion |
| `docs/known_issues.md` | Eight defects found during development (F1–F8): what was wrong, current fix status |
| `docs/pathology_score.md` | The planned (not yet built) morphometric pathology score |
| `CHANGELOG.md` | What changed, phase by phase |
| `config.example.yaml` | `AnalysisConfig` structure reference (not yet wired to `--config`) |

---

## Installation

> Requires **Python 3.11**. Conda is recommended.

```bash
conda create -n nerve_env python=3.11
conda activate nerve_env
pip install -r requirements.txt
```

> **tkinter** (GUI only) ships with Python but may need a separate step:
> - macOS Homebrew: `brew install python-tk`
> - conda: `conda install tk`

---

## Input files

Each analysis needs a **pair** of TIFF (or PNG/GIF for the mask) files:

| File | Description |
|------|-------------|
| TEM image | Raw greyscale electron micrograph |
| Mask image | Segmentation mask with integer pixel labels |

**Default mask pixel values** (configurable via `--myelin-val`/`--axoplasm-val`/`--mito-val`):

| Value | Tissue |
|-------|--------|
| 64 | Myelin sheath |
| 128 | Mitochondria |
| 192 | Axoplasm |
| 0 | Background |

**Supported file naming**:

| Pattern | TEM file | Mask file |
|---------|----------|-----------|
| A — underscore prefix | `tem_001.tif` or `axon_001.tif` | `mask_001.tif` |
| B — numeric prefix | `152. Axon 20K.tif` | `152. Mask 20K.tif` |

Pairs are matched by the leading integer (Pattern B) or the ID after the first `_`
(Pattern A). `--folder` (no `--recursive`) scans one directory non-recursively,
`.tif`/`.tiff` only. `--input-root`/`--groups`/`--recursive` scan recursively with an
extended glob (adds `.png`/`.gif`) and report unpaired files instead of silently
dropping them — use these for multi-subfolder datasets like this repo's own
`normal_data/`.

---

## Graphical interface

```bash
conda activate nerve_env
python measure_nerve_gui.py
```

A three-step workflow, left-hand rail:

**① Setup** — four sub-tabs, all scrollable:
- **Data** — one card per group (`+ Add group` for e.g. a future `treated` condition),
  each with a folder browser and a live pair-count preview.
- **Calibration** — pixel size (µm/px), with a derived `1 µm ≈ N px` readout, and a
  "pixel-space only" opt-out for uncalibrated data (see `CHANGELOG.md`'s "Not yet
  implemented" section for its current limits).
- **Segmentation** — a labeled **Auto-detect / Control / Pathology** mode selector
  (see `docs/known_issues.md` F3 for why this matters), mask-availability scanning,
  seven informational metric cards, and an advanced parameter accordion.
- **QC** — every threshold field, plus a live preview of flag counts recomputed
  against your last run as you edit (`docs/qc.md`).

Run is pinned below the sub-tabs — reachable regardless of which one is open.

**② Review** — two tabs:
- **Dashboard** — summary cards (images, axons, % flagged, mean g-ratio, mean
  mito occupancy) and a group-comparison boxplot+jitter strip.
- **Inspect** — TEM+mask overlay per image, with layer toggles, QC-status boundary
  coloring, a mitochondria sub-layer, click-to-inspect for any axon's full metric
  row, and Previous/Next navigation between images.

**③ Export** — output directory, a literal checklist of what will be written (with
row-count estimates and overwrite warnings), and an "Include 'mode' column" toggle
(off by default — see `docs/known_issues.md` F3).

The pre-redesign single-window GUI is retained as `measure_nerve_gui_legacy.py` if
you need it.

---

## Command-line interface

```
python measure_nerve.py [INPUT] [SCALE] [OPTIONS]
```

### Quick examples

```bash
# Single image — pixel size from scale bar, show plot
python measure_nerve.py \
    --tem "152. Axon 20K.tif" --mask "152. Mask 20K.tif" \
    --bar 191 1 --plot

# One folder, explicit pixel size
python measure_nerve.py --folder ./normal_data/152-159 --pixel-um 0.00524

# Multiple labeled groups in one run, with QC and mitochondrion-level output
python measure_nerve.py \
    --input-root . --groups normal_data:normal pathological_data:pathological \
    --pixel-um 0.00524 \
    --exclude-qc-failed --write-qc-report --write-mito-csv --write-group-csv
```

### Input (use one of these three)

| Option | Description |
|--------|-------------|
| `--folder DIR` | One folder, non-recursive, `.tif`/`.tiff` only |
| `--tem FILE --mask FILE` | A single image pair |
| `--input-root DIR --groups FOLDER:LABEL [FOLDER:LABEL ...]` | Multiple labeled groups, each scanned recursively |

`--group LABEL` (with `--folder` or `--tem`/`--mask`) stamps an explicit group label.
`--recursive` (with `--folder`) enables the recursive/extended-glob scan without
assigning a group.

### Scale (exactly one required)

| Option | Description |
|--------|-------------|
| `--pixel-um FLOAT` | Pixel size in µm/px |
| `--bar PX UM` | Scale bar length in pixels, then in µm |

### Segmentation

| Option | Default | Description |
|--------|---------|-------------|
| `--mode` | `auto` | `auto` / `normal` / `pathological` — see `docs/known_issues.md` F3 |
| `--myelin-threshold INT` | `200` | Myelin px count threshold for `auto` mode |
| `--myelin-val`, `--axoplasm-val`, `--mito-val INT` | `64`/`192`/`128` | Mask label values |
| `--smoothing-radius INT` | `1` | Morphological smoothing radius (px) |
| `--min-axon-area INT` | `200` | Minimum axon size to keep (px) |
| `--min-myelin-area INT` | `300` | Min myelin area for myelin metrics (px) |
| `--watershed-mode` | `weighted` | `weighted` / `simple` |
| `--watershed-weight` | `radius` | `radius` / `area` |
| `--watershed-compactness FLOAT` | `0.001` | Watershed compactness |
| `--watershed-beta FLOAT` | `1.0` | Size-bias strength for weighted watershed |
| `--assign-detached-myelin` | `none` | `none` / `nearest` |
| `--mito-hole-handling` | `fill` | `fill` (correct, default) / `legacy` (reproduces F1's bug — see `docs/known_issues.md`) |
| `--mito-assignment` | `centroid` | `centroid` (default) / `overlap` / `legacy` (reproduces F8's bug) |
| `--demyelination-reference FRACTION` | none | Reference `myelin_area_fraction_of_fov` for `image_demyelination_index`. No default — see `docs/known_issues.md` F4 |

### Quality control (see `docs/qc.md`)

| Option | Default | Description |
|--------|---------|-------------|
| `--exclude-qc-failed` | off | Mark (never drop) QC-failed axons |
| `--g-ratio-min`, `--g-ratio-max FLOAT` | `0.3`, `0.95` | |
| `--axon-circularity-min FLOAT` | `0.4` | |
| `--min-myelin-area-um2`, `--min-axon-area-um2`, `--min-fiber-area-um2 FLOAT` | none | `None` = not calibrated, flag stays off |
| `--min-mito-count-for-clustering INT` | `3` | |
| `--write-qc-report` | off | Also write `qc_report.csv` |

### Output

| Option | Default | Description |
|--------|---------|-------------|
| `--output-dir DIR` | input folder / `--input-root` | Where CSVs, plots, and the manifest go |
| `--plot` | off | Save + display overlay PNGs |
| `--write-mito-csv` | off | Also write `mitochondria_metrics.csv` (needs `--mito-assignment centroid` or `overlap`) |
| `--write-group-csv` | off | Also write `group_metrics.csv` |
| `--schema` | `full` | `full` (every column) / `legacy` (original pre-refactor column set + order) |

### Output files

| File | Contents |
|------|----------|
| `axons.csv` | One row per axon. Full column reference: `docs/metrics.md` §1–5 |
| `image_summary.csv` | One row per image. §7 |
| `mitochondria_metrics.csv` | One row per real mitochondrion (`--write-mito-csv`). §6 |
| `qc_report.csv` | Slim per-axon QC view (`--write-qc-report`). `docs/qc.md` |
| `group_metrics.csv` | One row per group (`--write-group-csv`). §8 |
| `run_manifest.json` | Always written: config, package versions, timestamp, inputs, skipped files, outputs, row counts, exclusion counts |
| `overlay_<id>.png` | Tissue overlay with axon IDs (`--plot`) |

---

## Computing pixel size

```
pixel_size_um = scale_bar_length_um / scale_bar_length_px
```

Example: a **1 µm** scale bar spanning **191 pixels** → `1 / 191 ≈ 0.00524 µm/px`.
Use ImageJ/Fiji to measure the scale bar length in pixels if it isn't already known.

---

## Troubleshooting

**No axons detected**
- Check that mask label values match the parameters (`--myelin-val`, `--axoplasm-val`, `--mito-val`)
- Try lowering `--min-axon-area`

**No TEM+mask pairs found**
- File names must contain `axon`, `tem`, or `mask` (case-insensitive), or start with a leading integer
- With `--input-root`/`--groups`, check `run_manifest.json`'s `skipped_files` for why a specific file wasn't paired

**Myelin metrics all `NaN` in normal tissue**
- Increase `--myelin-threshold` or set `--mode normal` explicitly
- Try lowering `--min-myelin-area`

**Myelin metrics all `NaN` in pathological tissue, and that seems expected but you want a demyelination number anyway**
- Supply `--demyelination-reference` — see `docs/known_issues.md` F4

**A mitochondrial number looks implausibly small, or `mito_outside_frac` is exactly `1.0`**
- Check `schema_version` in the row — `"1.0"` means `--mito-hole-handling legacy` was
  used (reproduces a known bug on purpose). Re-run with the default `fill`

**Too many spurious small fibers**
- Increase `--min-axon-area`

**Watershed boundaries look wrong**
- Try `--watershed-mode simple`
- Lower `--watershed-beta` for more equal boundary placement
