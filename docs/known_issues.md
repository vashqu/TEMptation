# Known issues and their fix status

Eight defects (F1–F8) were found auditing the pre-refactor tool (recorded in detail
in `AUDIT.md`, an internal engineering log). This document is the user-facing summary:
what was wrong, what it means for numbers you may already have, and whether it's
fixed today or still a documented limitation.

**If you have exported CSVs from before this tool had a `schema_version` column,
treat every mitochondrial number in them as unreliable and re-run.** See F1.

## F1 — mitochondrial area was ~2% of the true value (**fixed**, opt-in to revert)

The original mitochondrial-area calculation reached only ~1 pixel into each
mitochondrion, capturing a thin outer rim instead of the whole object. On a real
test image, this made reported mitochondrial area **2.1%** of the true value.
Every mitochondrial metric derived from area was affected: `mito_area_um2`,
`mito_density_per_um2`, `mito_mean_circularity`, `mito_mean_form_factor`,
`mito_mean_feret_um`, `mito_std_area_um2`, `mito_max_area_um2`, `mito_cv_area`,
`mito_area_skewness`, and (more subtly) `axon_area_um2` (biased low), `g_ratio`
(biased high), and every shape metric whose perimeter trace picked up the rim
boundary instead of real tissue edges.

**Fixed** as of `--mito-hole-handling fill`, the **default** since this fix landed.
The correct test is connectivity, not enclosure: a mitochondrion belongs to an axon
if it touches real axoplasm anywhere along its boundary, not only if it's fully
topologically enclosed (an earlier attempted fix using "fill enclosed holes" only
recovered 16.5% of true area — barely better than the original bug — because real
mitochondria often also touch myelin at one edge, disqualifying a fill-holes
approach entirely). `--mito-hole-handling legacy` reproduces the original bug
exactly, kept only for reproducing pre-fix numbers.

**`schema_version` in every output row tells you which algorithm produced it:**
`"2.0"` = fixed (`fill`), `"1.0"` = the original bug (`legacy`). Any pre-existing
analysis without this column predates the fix and its mitochondrial numbers should
be treated as unreliable.

## F2 — `mito_outside_frac` was structurally always `1.0` (**fixed**, same flag as F1)

`mito_outside_frac` (in `image_summary.csv`) was computed by comparing mitochondrial
pixels against the *raw* axoplasm mask, which by construction never overlaps
mitochondrial pixels — so the comparison found "outside" every single time, for
every image, regardless of actual segmentation quality. Both `mito_outside_area_um2`
and `mito_outside_frac` carried zero information under the original code.

**Fixed** under `--mito-hole-handling fill` (the default), which compares against
the corrected, mito-inclusive axoplasm mask instead. Still structurally `1.0` under
`--mito-hole-handling legacy`, by design — that mode exists only to reproduce old
numbers exactly.

## F3 — `mode` is a myelin-presence detector, not a group label (**documented, not a bug**)

`--mode auto` (the default) classifies each image independently by its own myelin
pixel count against `--myelin-threshold` (default 200px). This is **completely
independent** of the `group` label you assign per folder. An image inside your
`pathological_data` folder can still resolve to `mode=normal` if it happens to
contain enough myelin pixels — this was confirmed live: 4 of 6 images in a real
pathological dataset resolved as `mode=normal`.

This is not a bug — `mode` genuinely means "did this image's segmentation algorithm
find myelin," which is a real per-image fact, not a labeling error. But it's easy to
confuse with `group` given the shared vocabulary. The GUI's Segmentation tab
(Setup → Segmentation → "Segmentation mode") makes this distinction explicit with a
clearly-labeled selector, and the `mode` column is optional in GUI exports (off by
default — see Export → "Include 'mode' column"). See `docs/metrics.md` §1 for the
column's precise definition.

## F4 — myelin is unmeasurable in the pathological group by construction (**documented, not a bug**)

When myelin has fully detached from the axon (the typical pathological presentation),
`axon_vol_fraction` is identically `1.0` for every affected axon and
`myelin_area_um2`/`g_ratio`/`myelin_vol_fraction`/`mito_per_myelin` are all `NaN` —
verified on a real pathological dataset (25/25 axons). Both of CLAUDE.md's original
`demyelination_index` formula candidates reduce to `axon_area/fiber_area`, which is
blind precisely where it most needs to discriminate.

**Resolved at the image level**: `image_demyelination_index` is defined from
`myelin_area_fraction_of_fov` (a whole-mask quantity that survives detachment)
relative to a reference value — but that reference is a genuine biological
calibration choice this tool cannot derive on its own, so the column is `NaN`
unless you explicitly supply `--demyelination-reference` (see `docs/metrics.md`
§7.4). The per-axon `NaN`s themselves are expected and correctly tracked via the
non-excludable `qc_no_myelin` flag (`docs/qc.md`) — don't try to "fix" them by
excluding those axons; that would delete the pathology this tool exists to measure.

## F5 — `circularity` is not `1.0` for a perfect circle (**documented, unbiased alternative added**)

`skimage.measure.perimeter` overestimates a digitized boundary's true length, so the
legacy `circularity = 4π·area/perimeter²` column caps around 0.90–0.98 even for a
true disk, and is **size-dependent**: measured at 0.977 for a 10px-radius disk and
0.905 for a 160px-radius disk — the same shape, different apparent "circularity"
purely from raster size. Any README or documentation claiming "1 = perfect circle"
for this column was simply wrong; this has been corrected (`docs/metrics.md` §2).

The legacy column is kept exactly as-is (existing QC thresholds calibrated against
it remain valid) rather than "fixed" underneath users. `axon_circularity_crofton`
(and `axon_shape_irregularity_crofton`) use a materially less biased perimeter
estimator (`skimage`'s Crofton perimeter, ≈1.00 across the same size range where
legacy circularity drifts from 0.98 to 0.90) and are added **alongside**, not in
place of, the legacy columns. Prefer the Crofton variants for any comparison across
axons of different sizes.

## F6 — non-canonical mask pixel values (**detected, surfaced as a QC flag**)

At least one real mask file in this codebase's own test data contains 203 distinct
pixel values against the canonical `{0, 64, 128, 192}` — almost certainly a mask
export or resaving artifact, not intentional labeling.

Surfaced via `image_noncanonical_mask_frac` and `qc_mask_noncanonical` in
`image_summary.csv` (`docs/metrics.md` §7.5) — not silently ignored, not
auto-corrected. If `qc_mask_noncanonical` is `True` for an image, inspect that
mask's export pipeline before trusting its numbers.

## F7 — batch folder scanning silently dropped non-`.tif`/`.tiff` files (**fixed**)

The original folder scanner only globbed `*.tif`/`*.tiff`, and wasn't recursive into
subfolders at all — a real file in this codebase's own test data (`429. Mask 20K.gif`)
was silently dropped with no warning.

**Fixed**: `--input-root`/`--groups`/`--recursive` (and the GUI's Data tab) use a
recursive scan with an extended glob (`.tif`, `.tiff`, `.png`, `.gif`), and any file
that still can't be paired is reported as a `skipped_files` entry in
`run_manifest.json` and, in the GUI, directly on the group's card — not silently
dropped. The original non-recursive `--folder` behavior (without `--recursive`) is
kept exactly as-is for backward compatibility with existing scripts.

## F8 — mitochondrion double-counting at watershed boundaries (**fixed**), plus schema stability (**fixed**)

Two related findings:

1. The original per-axon mitochondria calculation cropped and re-labeled mitochondria
   independently for each fiber, so a single mitochondrion straddling a watershed
   boundary between two axons was fragmented into two pieces and counted (in part)
   toward both axons.
   **Fixed** via `--mito-assignment centroid` (the default): mitochondria are labeled
   once, globally, then each whole mitochondrion is assigned to exactly one axon (by
   which axon contains its centroid). `--mito-assignment overlap` (max pixel overlap)
   is offered as an alternative; `--mito-assignment legacy` reproduces the original
   crop-based double-counting exactly, for regression comparison only.
   `n_mito_assigned`/`n_mito_unassigned` (`image_summary.csv`) report how many
   mitochondria had no axon contact at all under the corrected modes.

2. An image with zero or one detected axons could previously produce a CSV missing
   columns present in other rows (e.g. `nearest_neighbor_um`, which needs ≥2 axons
   to be defined) — a schema that silently depended on how many objects happened to
   be found. **Fixed**: every output is reindexed to a fixed, versioned column set
   (`--schema {legacy,full}`) regardless of row count — a 0-axon image still emits
   the complete header, just with 0 rows.

## Column layout vs. value correctness — two independent axes

Two different flags control two different things, and they're easy to conflate:

- **`schema_version`** (a per-row CSV column, not a CLI flag) — which mitochondrial-
  area algorithm actually produced that row's *numbers* (F1 above). `"1.0"` = the
  original bug, `"2.0"` = fixed.
- **`--schema {legacy,full}`** (a CLI flag) — which *columns* get written at all.
  `legacy` writes only the original pre-refactor column set, in its original order,
  for downstream scripts that don't expect new columns. `full` (default) writes
  every column this version of the tool computes.

You can have `--schema legacy` with `schema_version="2.0"` (old columns, corrected
values) or `--schema full` with `schema_version="1.0"` (all columns, but
`--mito-hole-handling legacy` was used, so mitochondrial values reproduce the F1
bug on purpose, for regression testing). Check `schema_version` — not `--schema` —
when deciding whether a mitochondrial number is trustworthy.
