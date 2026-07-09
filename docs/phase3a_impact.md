# Phase 3a impact report — global mito-to-axon assignment (F8)

**Status: implemented and merged.** `--mito-assignment centroid` is now the CLI default (previous: `legacy`). `--mito-assignment legacy` remains available and is byte-exact against the Phase 0 golden baseline when combined with `--mito-hole-handling legacy` (verified).

## What was wrong

Per-axon mitochondrial metrics were computed by cropping to each watershed fiber's bounding box and intersecting with the raw mitochondria mask: `mito_crop & axon_mask_local`. A mitochondrion whose pixels straddle the watershed boundary between two adjacent axons would be split by this intersection — each axon's loop iteration would see only its own fragment as a separate connected component. The result: the same physical organelle counted as two mitochondria across two axons, its area split (not doubled) between them, and its shape metrics (circularity, form factor, Feret diameter) computed on a fragment rather than the true shape.

## What the fix does

`segmentation.assign_mito_to_axons` labels every mitochondrion in the (already Phase-2b-corrected) `mito_in_axon` mask **once, globally**, then assigns each whole organelle to exactly one axon: by the watershed label at its centroid pixel (default), falling back to maximum pixel overlap if the centroid lands outside every fiber (e.g. a concave mitochondrion whose centroid falls in myelin between its own arms — verified with a synthetic C-shaped test case). `metrics_mito.mito_metrics_from_regions` then computes the same aggregate metrics as before, but from these whole, correctly-assigned regions instead of re-labeled per-fiber fragments.

## Measured impact: none, on this dataset

Checked exhaustively before implementation: **0 of 1005 mitochondria across all 99 masks** in `normal_data`/`pathological_data` straddle a watershed boundary. Because of this, the crop-based and global-assignment algorithms are mathematically equivalent on every image in this dataset — provable directly: absent straddling, `mito_crop & axon_mask_local` for fiber X reduces to exactly `mito_in_axon ∩ (labels_ws == X)`, which is also what centroid/overlap assignment converges to when 100% of a mitochondrion's pixels already lie in one fiber's territory.

This was verified empirically, not just argued: `mito_metrics_for_axon` (legacy) and `mito_metrics_from_regions` (new) produce identical output (to 1e-9) for every axon in every image checked, including the full `normal_data/162-165` and `pathological_data` regression-anchor folders. `tests/golden_v2/` axon-level values are therefore **unchanged** by this phase. The only change to `tests/golden_v2/` is two new image-level diagnostic columns:

| column | meaning | value on this dataset |
|---|---|---|
| `n_mito_assigned` | mitochondria successfully attributed to an axon | matches total mito count per image |
| `n_mito_unassigned` | mitochondria in `mito_in_axon` with no fiber assignment (centroid *and* overlap both land outside every watershed region — e.g. attached only to axon_only tissue below the `min_axon_area_px` threshold, which never became a seeded fiber) | 0 for every image checked |

Both are `NaN` when `--mito-assignment legacy` is used (no global assignment is computed in that mode).

Confirmed on a full live run across the entire dataset (`--groups normal_data:normal pathological_data:pathological`, 100 images, 1165 axon rows): `n_mito_unassigned == 0` for all 100 images, 1012 mitochondria assigned, 0 unassigned.

## Why the fix ships anyway despite zero measured impact here

This is infrastructure correctness, not a here-and-now bug fix on this specific dataset: it removes an algorithmic edge case (straddling mitochondria) that this dataset happens not to exercise, but that other/future datasets — different magnification, different segmentation model, more crowded axon fields — plausibly could. It is also the foundation Phase 3b's per-mitochondrion CSV (`mitochondria_metrics.csv`) needs: a table with "one row per real mitochondrion" is only meaningful once no mitochondrion can appear as two rows.

## Limitations

1. Verified only on this dataset (99 masks, all `20K` magnification, all TEM axon micrographs from the same acquisition pipeline). The straddling scenario this fix addresses is untested on real straddling data — only on synthetic cases (see `tests/test_segmentation_mito_assignment.py`).
2. `n_mito_unassigned` being 0 everywhere in this dataset does not mean it will always be 0. It is a genuine QC signal and should be watched on new data.
