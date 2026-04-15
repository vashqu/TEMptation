# Nerve Fiber Analyzer

Segment and measure nerve fibers from **TEM images + segmentation masks**.  
Extracts axon morphology, myelin metrics, mitochondria shape features, and image-level summaries.

Two interfaces are provided — a **GUI** and a **CLI** — sharing the same core engine.

---

## Features

- **Axon-level metrics:** area, diameters, g-ratio, myelin thickness, circularity, convexity, AVF/MVF
- **Mitochondria features:** count, area, density, circularity, form factor, Feret diameter, spatial clustering
- **Image-level summary:** fiber density, mean g-ratio, myelin fraction, mitochondria outside axons
- **Two watershed modes:** simple or size-weighted (larger axons claim proportionally more myelin territory)
- **Auto mode detection:** automatically identifies normal vs. pathological tissue from myelin pixel count
- **Batch processing:** process an entire folder of image pairs in one run (GUI and CLI)

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

## Input Files

Each analysis requires a **pair** of TIFF files:

| File | Description |
|------|-------------|
| TEM image | Raw greyscale electron micrograph |
| Mask image | Segmentation mask with integer pixel labels |

**Default mask pixel values** (configurable):

| Value | Tissue |
|-------|--------|
| 64 | Myelin sheath |
| 128 | Mitochondria |
| 192 | Axoplasm |
| 0 | Background |

**Supported file naming** for batch/folder mode:

| Pattern | TEM file | Mask file |
|---------|----------|-----------|
| A — underscore prefix | `tem_001.tif` or `axon_001.tif` | `mask_001.tif` |
| B — numeric prefix | `152. Axon 20K.tif` | `152. Mask 20K.tif` |

Pairs are matched by the leading integer (Pattern B) or the ID after the first `_` (Pattern A).

---

## Graphical Interface

```bash
conda activate nerve_env
python measure_nerve_gui.py
```

### Single image mode (default)

1. Leave **Input mode** set to **Single image**
2. Browse for your **TEM image** (optional — needed only for the overlay plot)
3. Browse for your **Mask image** (required)
4. Adjust **Pixel size (µm/px)** — see [Computing pixel size](#computing-pixel-size)
5. Click **▶ Run Analysis**
6. View results in the **Axon Metrics** and **Image Summary** tabs
7. Click **Show Overlay Plot** to visualise the segmentation
8. Click **Export CSV…** to save results

### Batch folder mode

1. Select **Batch folder** in the Input mode row
2. Click **Browse folder…** and select the folder with your image pairs
3. The info line shows how many pairs were detected
4. Set parameters and click **▶ Run Analysis**
5. The result tables show combined data from all images (`image_id` column identifies each)
6. Click **Export CSV…** to save the aggregated results

---

## Command-Line Interface

```
python measure_nerve.py [INPUT] [SCALE] [OPTIONS]
```

### Quick examples

```bash
# Single image — pixel size from scale bar, show plot
python measure_nerve.py \
    --tem "152. Axon 20K.tif" \
    --mask "152. Mask 20K.tif" \
    --bar 191 1 \
    --plot

# Whole folder — explicit pixel size
python measure_nerve.py \
    --folder ./normal_data/152-159 \
    --pixel-um 0.00524

# Pathological tissue — simple watershed, save to custom folder
python measure_nerve.py \
    --tem img.tif --mask mask.tif \
    --bar 191 1 \
    --mode pathological \
    --watershed-mode simple \
    --output-dir ./results
```

### All options

| Option | Default | Description |
|--------|---------|-------------|
| `--folder DIR` | — | Folder to process (batch mode) |
| `--tem FILE` | — | TEM image (single mode) |
| `--mask FILE` | — | Mask image (single mode) |
| `--pixel-um FLOAT` | — | Pixel size in µm/px |
| `--bar PX UM` | — | Scale bar length in pixels then µm |
| `--mode` | `auto` | `auto` / `normal` / `pathological` |
| `--myelin-threshold INT` | `200` | Myelin px count threshold for auto mode |
| `--myelin-val INT` | `64` | Mask value for myelin |
| `--axoplasm-val INT` | `192` | Mask value for axoplasm |
| `--mito-val INT` | `128` | Mask value for mitochondria |
| `--smoothing-radius INT` | `1` | Morphological smoothing radius (px) |
| `--min-axon-area INT` | `200` | Minimum axon size to keep (px) |
| `--min-myelin-area INT` | `300` | Min myelin area for myelin metrics (px) |
| `--watershed-mode` | `weighted` | `weighted` / `simple` |
| `--watershed-weight` | `radius` | `radius` / `area` |
| `--watershed-compactness FLOAT` | `0.001` | Watershed compactness |
| `--watershed-beta FLOAT` | `1.0` | Size-bias strength for weighted watershed |
| `--assign-detached-myelin` | `none` | `none` / `nearest` |
| `--output-dir DIR` | input folder | Where to save CSVs and plots |
| `--plot` | off | Save and show overlay plots |

### Output files

| File | Contents |
|------|----------|
| `axons.csv` | One row per axon per image |
| `image_summary.csv` | One row per image |
| `overlay_<id>.png` | Tissue overlay with axon IDs (when `--plot`) |

---

## Computing pixel size

```
pixel_size_um = scale_bar_length_um / scale_bar_length_px
```

Example: a **1 µm** scale bar spanning **191 pixels** → `1 / 191 ≈ 0.00524 µm/px`

Use ImageJ/Fiji to measure the scale bar length in pixels if it is not already known.

---

## Output columns (axons.csv)

| Column | Description |
|--------|-------------|
| `image_id` | Source image identifier |
| `mode` | `normal` or `pathological` |
| `axon_id` | Unique axon label within the image |
| `axon_area_um2` | Axon cross-sectional area (µm²) |
| `myelin_area_um2` | Myelin area (µm²); NaN in pathological |
| `fiber_area_um2` | Total fiber area (axon + myelin) |
| `d_inner_um` | Equivalent axon diameter (µm) |
| `d_outer_um` | Equivalent fiber diameter (µm); NaN in pathological |
| `g_ratio` | d_inner / d_outer; NaN in pathological |
| `myelin_thickness_um` | (d_outer − d_inner) / 2; NaN in pathological |
| `perimeter_um` | Axon perimeter (µm) |
| `eccentricity` | 0 = circle, 1 = line |
| `solidity` | area / convex_hull_area |
| `circularity` | 4π·area / perimeter² (1 = perfect circle) |
| `convexity` | convex_hull_perimeter / perimeter |
| `axon_vol_fraction` | Axon area / fiber area (AVF) |
| `myelin_vol_fraction` | Myelin area / fiber area (MVF); NaN in pathological |
| `centroid_x_um`, `centroid_y_um` | Centroid in µm |
| `centroid_x_px`, `centroid_y_px` | Centroid in pixels |
| `mito_count` | Number of mitochondria |
| `mito_area_um2` | Total mitochondria area |
| `mito_density_per_um2` | Mitochondria count per µm² |
| `mito_mean_circularity` | Mean circularity of individual mito objects |
| `mito_mean_form_factor` | Mean form factor (perimeter²/4π·area) |
| `mito_mean_feret_um` | Mean max Feret diameter (µm) |
| `mito_std_area_um2` | Std of mito areas (size heterogeneity) |
| `mito_max_area_um2` | Largest single mito area |
| `mito_cv_area` | Coefficient of variation of mito areas |
| `mito_area_skewness` | Skewness of mito area distribution |
| `mito_mean_dist_centroid_um` | Mean distance of mito centroids from axon centroid |
| `mito_std_dist_centroid_um` | Std of those distances |
| `nearest_neighbor_um` | Distance to nearest axon centroid (µm) |

---

## Troubleshooting

**No axons detected**
- Check that mask label values match the parameters (`--myelin-val`, `--axoplasm-val`, `--mito-val`)
- The GUI displays detected pixel values when you load a mask
- Try lowering `--min-axon-area`

**No TEM+mask pairs found (batch mode)**
- File names must contain `axon`, `tem`, or `mask` (case-insensitive)
- Pairs are matched by leading integer or by the ID after the first `_`

**Myelin metrics all NaN in normal tissue**
- Increase `--myelin-threshold` or set `--mode normal` explicitly
- Try lowering `--min-myelin-area`

**Too many spurious small fibers**
- Increase `--min-axon-area`

**Watershed boundaries look wrong**
- Try `--watershed-mode simple`
- Lower `--watershed-beta` for more equal boundary placement
