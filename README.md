# TEMptation

TEMptation measures nerve fibers in transmission electron microscopy (TEM)
images from a raw TEM image and a pixel-labeled segmentation mask. It can also
optionally call AxonDeepSeg to create an automatic mask from an unlabeled TEM
image before running the same morphometry pipeline.

Two interfaces are available:

- `measure_nerve.py` - command-line interface for scripts and batch runs
- `measure_nerve_gui.py` - Tkinter GUI for interactive use

## What It Measures

Axon-level outputs include axon area, fiber area, equivalent diameters,
g-ratio, myelin thickness, perimeter, eccentricity, solidity, circularity,
convexity, nearest-neighbor distance, and mitochondria features when
mitochondria are labeled.

Image-level outputs include fiber count, field-of-view area, fiber density,
mean g-ratio, mean myelin thickness, myelin area fraction, and mitochondria
outside axon territory.

Note: the axon-level `axon_vol_fraction` and `myelin_vol_fraction` columns are
per-fiber area fractions computed from each fiber's axon and myelin area. They
are preserved for compatibility with existing outputs.

## Installation

The lightweight environment is enough when you already have TEM images and
pixel-labeled masks.

```bash
conda create -n nerve_env python=3.11
conda activate nerve_env
pip install -r requirements.txt
```

Verify the install:

```bash
python -c "import numpy, scipy, skimage, pandas, tifffile, matplotlib; print('OK')"
```

For the GUI, Tkinter may need a separate install depending on your Python:

```bash
# conda
conda install tk

# macOS Homebrew Python
brew install python-tk
```

### Optional AxonDeepSeg Environment

Automatic segmentation needs AxonDeepSeg and its deep-learning dependencies.
Create the optional environment from the TEMptation directory:

```bash
conda env create -f environment-ads.yml
conda activate temptation_ads
```

You can also keep AxonDeepSeg in a separate environment and point TEMptation to
that Python:

```bash
export AXONDEEPSEG_PYTHON=/path/to/axondeepseg/env/bin/python
```

TEMptation automatically looks for a sibling `../axondeepseg` clone. Override
that with `AXONDEEPSEG_PATH` or the CLI option `--ads-package-dir`.

## Input Masks

The standard workflow expects a raw TEM image and a label mask with the same
height and width.

Default mask values:

| Value | Meaning |
| --- | --- |
| `0` | Background |
| `64` | Myelin |
| `128` | Mitochondria |
| `192` | Axoplasm |

The label values are configurable in both the CLI and GUI.

Batch folder mode supports these naming patterns:

| TEM image | Mask image | Matched ID |
| --- | --- | --- |
| `tem_001.tif` or `axon_001.tif` | `mask_001.tif` | `001` |
| `152. Axon 20K.tif` | `152. Mask 20K.tif` | `152` |

## CLI Usage

Run morphometry from an existing mask:

```bash
python measure_nerve.py \
  --tem "152. Axon 20K.tif" \
  --mask "152. Mask 20K.tif" \
  --pixel-um 0.00524 \
  --output-dir ./results
```

Use a scale bar instead of an explicit pixel size:

```bash
python measure_nerve.py \
  --tem "152. Axon 20K.tif" \
  --mask "152. Mask 20K.tif" \
  --bar 191 1
```

Process a folder of matched TEM/mask pairs:

```bash
python measure_nerve.py \
  --folder ./normal_data/152-159 \
  --pixel-um 0.00524 \
  --output-dir ./results
```

Run AxonDeepSeg on an unlabeled TEM image and then measure it:

```bash
python measure_nerve.py \
  --tem "152. Axon 20K.tif" \
  --auto-segment \
  --ads-model generalist \
  --pixel-um 0.00524 \
  --output-dir ./results
```

Save only the TEM image and generated mask for manual refinement:

```bash
python measure_nerve.py \
  --tem "152. Axon 20K.tif" \
  --auto-segment \
  --auto-segment-action save \
  --ads-model unmyelinated-TEM \
  --output-dir ./manual_refinement
```

Useful options:

| Option | Description |
| --- | --- |
| `--mode auto|normal|pathological` | Choose myelin-aware or axon-only metrics; `auto` uses myelin pixel count. |
| `--watershed-mode weighted|simple` | Assign myelin territory with size-weighted or simple watershed. |
| `--assign-detached-myelin nearest` | Attribute detached myelin pixels to nearest axon for myelin metrics. |
| `--plot` | Save an overlay PNG with tissue labels and axon IDs. |
| `--ads-model generalist|unmyelinated-TEM` | Choose the AxonDeepSeg model for automatic segmentation. |
| `--ads-python PATH` | Use a specific Python executable that has AxonDeepSeg installed. |
| `--ads-model-path PATH` | Use an already downloaded AxonDeepSeg model folder. |

Run `python measure_nerve.py --help` for the full option list.

## GUI Usage

Start the GUI:

```bash
python measure_nerve_gui.py
```

Single image mode:

1. Select `Single image`.
2. Choose the TEM image and mask image.
3. Set `Pixel size (µm/px)`.
4. Click `Run Analysis`.
5. Review the `Axon Metrics` and `Image Summary` tabs.
6. Export CSVs when needed.

Auto segment mode:

1. Select `Auto segment`.
2. Choose an unlabeled TEM image.
3. Choose `generalist` or `unmyelinated-TEM`.
4. Choose `analyze` to measure immediately or `save` to save the generated
   mask for manual refinement.
5. Click `Run Analysis`.

AxonDeepSeg automatic masks do not include mitochondria labels at this point.
Mitochondria metrics require manual refinement with label value `128`.

Batch folder mode:

1. Select `Batch folder`.
2. Choose a folder containing matched TEM/mask pairs.
3. Set the pixel size and parameters.
4. Click `Run Analysis`.
5. Export the combined CSV outputs.

## Outputs

The CLI writes outputs to `--output-dir` or to the input folder by default.

| File | Contents |
| --- | --- |
| `axons.csv` | One row per detected axon/fiber. |
| `image_summary.csv` | One row per image. |
| `overlay_<id>.png` | Optional overlay when `--plot` is used. |
| `*_ads_mask.tif` | TEMptation-compatible mask generated by AxonDeepSeg. |

## Troubleshooting

No axons detected:

- Confirm the mask values match the configured labels.
- Lower `--min-axon-area` if small axons are being filtered.
- Check that TEM and mask dimensions match.

Myelin metrics are `NaN`:

- The image may have been resolved as pathological mode.
- Use `--mode normal` if myelin should be measured.
- Lower `--min-myelin-area` if valid myelin regions are small.

AxonDeepSeg does not run:

- Confirm AxonDeepSeg is installed in the current environment or set
  `AXONDEEPSEG_PYTHON`.
- Use `--ads-model-path` if models were downloaded manually.
- CPU inference can be slow; use `--ads-gpu-id 0` when a compatible GPU is
  available.

Slow import or font-cache warnings:

- Use a writable Matplotlib cache directory, for example:
  `export MPLCONFIGDIR=/tmp/matplotlib-cache`.
- Non-plot analysis no longer imports `matplotlib.pyplot` at module import time.

## Development

Before pushing changes, run a quick syntax check and at least one representative
CLI smoke test with local sample data:

```bash
python -m py_compile auto_segment.py measure_nerve.py measure_nerve_gui.py
python measure_nerve.py --tem "../normal_data/152-159/152. Axon 20K.tif" --mask "../normal_data/152-159/152. Mask 20K.tif" --pixel-um 0.00524 --output-dir ./results_TEMPORARY
```

Generated CSVs, overlays, logs, `__pycache__/`, local virtual environments, and
`*_TEMPORARY/` result folders are ignored by `.gitignore`.

## Citation

If you use the automatic segmentation option, cite AxonDeepSeg:

Zaimi, A., Wabartha, M., Herman, V. et al. AxonDeepSeg: automatic axon and
myelin segmentation from microscopy data using convolutional neural networks.
Sci Rep 8, 3816 (2018). https://doi.org/10.1038/s41598-018-22181-4
