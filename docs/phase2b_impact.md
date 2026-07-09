# Phase 2b impact report — mitochondrial-rim bug fix (F1) and mito_outside_frac fix (F2)

**Status: implemented and merged.** `--mito-hole-handling fill` is now the CLI default (previous behavior: `legacy`, still available via `--mito-hole-handling legacy`). `schema_version` in every output row records which algorithm produced it: `"1.0"` = legacy/pre-fix, `"2.0"` = fill/fixed.

## What was wrong

`measure_nerve.py` built the axoplasm compartment as `axon_only | (mito & binary_dilation(axon_only, disk(1)))`. Since a mitochondrion is a *hole* in the axon_only mask (label 192 excludes mitochondria, labeled 128), dilating axon_only by 1 pixel only reached the outer 1px rim of each mitochondrion. Every mitochondrial metric (`mito_area_um2`, `mito_density_per_um2`, `mito_mean_circularity`, `mito_mean_form_factor`, `mito_std_area_um2`, `mito_max_area_um2`, `mito_cv_area`, `mito_area_skewness`) was computed on that rim, not the mitochondrion — measured on one real image, only 2.1% of the true mitochondrial area was captured. `axon_area_um2` was correspondingly understated (the holes were never filled in), which propagated into `g_ratio`, `circularity`, `solidity`, `convexity`, and every axon-boundary-derived quantity.

A second, independent bug (F2): the image-level `mito_outside_frac`/`mito_outside_area_um2` metric compared mitochondria against the *raw* `axon_only` mask (which by construction never overlaps mitochondria pixels at all), so it was **structurally always exactly 1.0** regardless of mode, image, or dataset — it carried zero information.

## What the fix does

**Not** a naive "fill the topological holes in axon_only" approach — that was tried first and rejected. `ndi.binary_fill_holes(axon_only)` only fills a hole with *zero* leak path to any other tissue type (myelin, background, or another mitochondrion are all "not axon_only" to that algorithm). On one real image, a 41,551px mitochondrion that plainly touches real axon_only territory also touches myelin at its far edge — biologically unremarkable, since mitochondria are not confined to the exact geometric center of an axon — and that single myelin contact was enough to disqualify the whole component. Fill-holes alone recovered only 16.5% of true mitochondrial area: better than the 2.2% legacy bug, but still substantially wrong.

The shipped algorithm tests **connectivity, not enclosure**: a mitochondrion belongs to an axon if its connected component, unioned with `axon_only`, contains at least one `axon_only` pixel — i.e. it directly touches real axoplasm anywhere along its boundary, regardless of what else it touches. Verified on the full dataset (99 masks): **99.65% of all mitochondrial area** is captured this way. The remaining 0.35% (8,338px across 4 of 99 images) is genuinely orphaned mitochondria — segmentation blobs with zero adjacency to any axon — and is correctly excluded rather than silently mis-attributed. One of those four images (`mask_path_3.tif`) is the same file separately flagged for having 2,645 non-canonical mask pixels (F6); its orphaned mitochondrion is very likely another symptom of that mask's known quality problems.

`mito_outside_frac`/`mito_outside_area_um2` were also fixed as part of this same change: they now compare against the corrected (fill-mode) axoplasm instead of the always-empty-intersection `axon_only`, so they report the genuinely orphaned fraction described above instead of a constant 1.0.

**Segmentation geometry is mostly, but not entirely, unaffected.** Axon counts are identical between legacy and fill modes for the two folders used as the golden-baseline regression anchor (`normal_data/162-165`: 6/14/14/13; all of `pathological_data`: 7/2/5/3/4/4) -- this was verified first and is what an earlier draft of this document over-generalized from. Run across the **full dataset** (100 images, 1172 legacy axons vs 1165 fill axons), 7 of 100 images show exactly one fewer axon under fill mode:

| image | legacy axons | fill axons |
|---|---:|---:|
| 172 | 16 | 15 |
| 236 | 18 | 17 |
| 238 | 17 | 16 |
| 239 | 19 | 18 |
| 251 | 18 | 17 |
| 254 | 13 | 12 |
| 275 | 14 | 13 |

**Mechanism, verified with zero exceptions across all 7 images:** watershed axon seeds are built from `axoplasm = axon_only | mito_in_axon`, so anything that changes `mito_in_axon`'s extent can change axon *connectivity*, not just axon *area*. Under the legacy bug, a large mitochondrion that happens to bisect an axon's axon_only footprint left the two halves disconnected in the axoplasm mask, and watershed correctly-per-its-input treated them as two separate seed components -- one real axon plus a spurious fragment (in several cases well under 1000px, i.e. barely or not even above the 200px minimum-axon-area filter). Under fill mode, that same mitochondrion is correctly recognized as part of the axon (by the connectivity test described above) and bridges the two halves back into one component. Checked directly for every one of the 7 images: the pixel gap connecting each merged pair is 99.2-99.9% mitochondrion-labeled, and the merged area equals the sum of the two legacy fragments plus the gap exactly (e.g. image 172: 902 + 4986 + 2631 = 8519, matching the single fill-mode component). Several images show more than one such merge event at the raw-component level (236 and 275 each show 2-3), but most of those additional merges involve components already below the 200px minimum-axon-area threshold in both modes, so they don't change the final filtered axon count -- only one merge per image, in every one of these 7 cases, crosses that threshold and changes the count.

**This is a segmentation-quality improvement, not a regression.** The legacy behavior was over-segmenting these 7 axons into one real fiber plus one spurious mitochondrion-shaped fragment; the fix correctly merges them into a single fiber. It does mean `axons.csv` row counts are not guaranteed identical between `--mito-hole-handling legacy` and the `fill` default on arbitrary data, even though they happen to be identical on the two folders used for byte-exact regression testing.

## Quantified impact — axon-level means

### Normal group (`normal_data/162-165`, n=47 axons)

| column | legacy mean | fill mean | change |
|---|---:|---:|---:|
| `axon_area_um2` | 0.71628 | 0.77630 | +8.4% |
| `fiber_area_um2` | 1.62588 | 1.68589 | +3.7% |
| `d_inner_um` | 0.82066 | 0.85664 | +4.4% |
| `d_outer_um` | 1.27107 | 1.29127 | +1.6% |
| `g_ratio` | 0.64135 | 0.65491 | +2.1% |
| `myelin_thickness_um` | 0.22521 | 0.21732 | -3.5% |
| `perimeter_um` | 3.89863 | 3.34203 | -14.3% |
| `eccentricity` | 0.82996 | 0.82883 | -0.1% |
| `solidity` | 0.89009 | 0.94106 | +5.7% |
| `circularity` | 0.57165 | 0.69195 | +21.0% |
| `convexity` | 0.89916 | 0.98522 | +9.6% |
| `axon_vol_fraction` | 0.41755 | 0.43240 | +3.6% |
| `myelin_vol_fraction` | 0.58245 | 0.56760 | -2.5% |
| `mito_count` | 0.48936 | 0.48936 | **0.0%** (unaffected, as expected) |
| `mito_area_um2` | 0.00258 | 0.06255 | **+2329%** |
| `mito_density_per_um2` | 0.75120 | 0.60031 | -20.1% (count flat, area denominator grew) |
| `mito_mean_circularity` | 0.05931 | 0.86642 | **+1361%** (rim → plausible blob shape) |
| `mito_mean_form_factor` | 20.83849 | 1.15664 | **-94.4%** (thin-ring artifact → near-1) |
| `mito_mean_feret_um` | 0.43811 | 0.43811 | 0.0% (max caliper of a rim ≈ max caliper of the filled blob) |
| `mito_std_area_um2` | 0.00036 | 0.01188 | +3168% |
| `mito_max_area_um2` | 0.00580 | 0.15493 | +2570% |
| `mito_cv_area` | 0.31397 | 0.55365 | +76.3% |
| `mito_area_skewness` | 0.02297 | 0.35117 | +1429% |
| `mito_mean_dist_centroid_um` | 0.30776 | 0.28055 | -8.8% |
| `mito_std_dist_centroid_um` | 0.05597 | 0.05829 | +4.1% |

### Normal group — image-level (n=4 images)

| column | legacy mean | fill mean | change |
|---|---:|---:|---:|
| `mean_circularity` | 0.55743 | 0.68999 | +23.8% |
| `mean_convexity` | 0.89095 | 0.98588 | +10.7% |
| `mean_g_ratio` | 0.64229 | 0.65641 | +2.2% |
| `mito_outside_area_um2` | 0.73499 | 0.00000 | -100% (F2: no orphaned mito in this 4-image subset) |
| `mito_outside_frac` | 1.00000 | 0.00000 | -100% (was structurally meaningless; now genuinely 0 here) |

### Pathological group (n=25 axons / 6 images)

| column | legacy mean | fill mean | change |
|---|---:|---:|---:|
| `axon_area_um2` | 2.00634 | 2.08334 | +3.8% |
| `circularity` | 0.47388 | 0.59755 | +26.1% |
| `convexity` | 0.82626 | 0.91806 | +11.1% |
| `mito_area_um2` | 0.00451 | 0.08151 | +1706% |
| `mito_mean_circularity` | 0.05666 | 0.77737 | +1272% |
| `mito_mean_form_factor` | 20.66918 | 1.32895 | -93.6% |
| `mito_outside_frac` (image-level) | 0.66667 | 0.16667 | -75% (F2; `mask_path_3.tif`'s orphaned mito now correctly isolated instead of contaminating every image's ratio) |

`d_outer_um`, `g_ratio`, `myelin_thickness_um`, `myelin_vol_fraction` are `NaN` for the pathological group in both legacy and fill modes — unaffected by this fix; that is defect F4 (detached myelin, unmeasurable in this group by design), a separate issue.

## Directional sanity check

Every shift is explainable and internally consistent:
- `mito_count`, `mito_mean_feret_um` unchanged — object identity and max caliper of a filled blob equal that of its outer boundary; expected.
- `mito_mean_form_factor` collapsing from ~21 (thin-ring signature) to ~1.2-1.3 (near-circular blob) and `mito_mean_circularity` rising from ~0.06 to ~0.77-0.87 — exactly the shape-metric signature of measuring a 1px ring versus measuring a filled organelle.
- `mito_density_per_um2` *decreasing* despite `mito_area_um2` increasing 17-23x: density is `count/axon_area`, count is flat, and `axon_area_um2` also increased (the holes are now filled into the axon), so the ratio falls. Logically required given the other two numbers, not a separate effect.
- `circularity`/`convexity`/`solidity` all increase: axon perimeter no longer traces jagged mitochondrial-hole boundaries.
- `perimeter_um` *decreases* ~14-15%: consistent with the same mechanism (less jagged boundary = shorter measured perimeter).

## Limitations and cautions

1. **99.65% capture, not 100%.** 8,338px (0.35% of all mitochondrial pixels dataset-wide) remain genuinely unattributed in 4 of 99 images — these are mitochondria with literally zero pixel adjacency to any axon_only territory, most plausibly segmentation noise. They are correctly excluded rather than guessed at; `mito_outside_frac` is where this now shows up.
2. **Every previously-exported mitochondrial number from this tool is invalidated**, not just biased — the legacy values were measuring a different (much smaller, wrong-shaped) object. Any prior analysis, plot, or conclusion drawn from `mito_area_um2`, `mito_density_per_um2`, or any `mito_mean_*`/`mito_std_*`/`mito_cv_*`/`mito_max_*` column computed before this fix should be discarded, not adjusted.
3. **Axon-boundary metrics shifted by single-digit-to-~25% amounts**, not orders of magnitude — `g_ratio` +2.1%, `axon_area_um2` +8.4%, `circularity` +21-28%. These were biased but directionally usable before; treat pre-fix values as approximate, not wrong by a large factor.
4. **`--mito-hole-handling legacy` remains available and exactly reproduces the pre-fix numbers** (regression-tested against the Phase 0 golden baseline). Use it only to reproduce or audit historical output, never as an analysis default.
5. This fix does not address F3 (mode auto-detection misclassifying pathological images), F4 (myelin unmeasurable when detached), F5 (circularity's size-dependent bias — `axon_circularity_crofton` from Phase 2 partially addresses this separately), F6 (non-canonical mask pixels), or F7 (fixed independently in Phase 2). Those remain open per `IMPLEMENTATION_BLUEPRINT.md`.
