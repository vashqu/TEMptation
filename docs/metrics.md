# Metric dictionary

Every column TEMptation can produce, generated from the code (`temptation/schema.py`,
`metrics_axon.py`, `metrics_mito.py`, `metrics_spatial.py`, `summaries.py`, `pipeline.py`),
not written from memory. If a column here doesn't match what you see in a CSV, the code
is authoritative — please file that as a docs bug.

**Units.** `_um`/`_um2` columns are already calibrated by `pixel_size_um`; `_px` columns
are raw pixels. Ratios (`g_ratio`, `circularity`, `mito_occupancy_ratio`, ...) are
dimensionless and identical in both unit systems.

**NaN conventions.** A metric is `NaN` when it is *undefined* for that row (e.g. g-ratio
for a pathological axon with no myelin, or any std/skewness statistic computed from too
few objects) — never as a stand-in for zero. `mito_count = 0` and `mito_occupancy_ratio = 0`
are both real, meaningful zeros; `mito_fragmentation_index` is `NaN` at zero mitochondria
because "how fragmented is nothing" is undefined, not zero.

**Not diagnostic.** Every metric below is a morphometric measurement derived from a
segmentation mask. None of it is validated for, or intended for, clinical diagnosis.
See `docs/pathology_score.md` for the fuller statement.

---

## 1. Identity & metadata columns

Present on every row of `axons.csv` and `image_summary.csv`, stamped by
`pipeline.analyze_dataset` (not computed per-metric).

| metric_name | level | formula | interpretation | notes |
|---|---|---|---|---|
| `image_id` | axon, image | — | Source image identifier (from filename) | |
| `axon_id` | axon | — | Watershed label of this axon within its image | Unique only within one image, not globally |
| `mode` | axon, image | `resolve_mode(mask)`: `"normal"` if myelin px count ≥ `--myelin-threshold`, else `"pathological"` | Which segmentation algorithm variant produced this row (normal computes g-ratio/myelin; pathological is axon-only) | **Not the same as `group`.** Auto-detected per image from mask content; two images in the same `group` can resolve to different `mode`. See `docs/known_issues.md` F3 |
| `group` | axon, image | user-supplied | Experimental condition label (`normal`, `pathological`, `treated`, ...) | Set via `--group` or `--groups FOLDER:LABEL`; independent of `mode` |
| `image_path`, `mask_path` | axon, image | — | Resolved source file paths | |
| `pixel_size_um` | axon, image | user-supplied | Calibration used for this run | See README "Computing pixel size" |
| `schema_version` | axon, image | `"2.0"` if `--mito-hole-handling fill` (default), else `"1.0"` | Which mitochondrial-area algorithm produced this row's numbers | **Not the same as `--schema`**, which controls column *layout* (legacy vs full), not value correctness. See `docs/known_issues.md` |

---

## 2. Core axon geometry & shape (`axons.csv`)

`fiber_area_um2 = axon_area_um2 + myelin_area_um2` under the default
`--assign-detached-myelin none`; the equivalent-diameter g-ratio convention is used throughout.

| metric_name | level | formula | interpretation | notes |
|---|---|---|---|---|
| `axon_area_um2` | axon | pixel count × `pixel_size_um²` | Axoplasm cross-sectional area | Mito-inclusive under `mito_hole_handling=fill` (default) — see `docs/known_issues.md` F1 |
| `myelin_area_um2` | axon | pixel count × `pixel_size_um²` | Myelin sheath area for this fiber | `NaN` in pathological mode or when below `--min-myelin-area` |
| `fiber_area_um2` | axon | watershed region area × `pixel_size_um²` | Total fiber area (axon + myelin) | |
| `d_inner_um` | axon | `2·√(axon_area_um2 / π)` | Equivalent axon diameter | |
| `d_outer_um` | axon | `2·√(fiber_area_um2 / π)` | Equivalent fiber diameter | `NaN` in pathological mode |
| `g_ratio` | axon | `d_inner_um / d_outer_um` | Axon-to-fiber diameter ratio; lower = more myelin relative to axon | `NaN` in pathological mode. Typical healthy range ≈ 0.6–0.8 |
| `myelin_thickness_um` | axon | `(d_outer_um − d_inner_um) / 2` | Equivalent myelin sheath thickness | `NaN` in pathological mode |
| `axon_vol_fraction` | axon | `axon_area_um2 / fiber_area_um2` | AVF — fraction of the fiber that is axoplasm | Identically `1.0` for every pathological axon (myelin unmeasurable — F4) |
| `myelin_vol_fraction` | axon | `myelin_area_um2 / fiber_area_um2` | MVF — fraction of the fiber that is myelin | `NaN` in pathological mode |
| `perimeter_um` | axon | `skimage.measure.perimeter` × `pixel_size_um` | Axon boundary length | Digital-boundary estimator; see `circularity` note below |
| `eccentricity` | axon | ellipse fit eccentricity | `0` = circle, →`1` = elongated line | |
| `solidity` | axon | `area / convex_hull_area` | `1` = fully convex; lower = concave/irregular boundary | |
| `convexity` | axon | `convex_hull_perimeter / perimeter` | `1` = convex; lower = irregular boundary | |
| `circularity` | axon | `4π·area / perimeter²` | Shape roundness | **Not `1.0` for a perfect circle** — `skimage.measure.perimeter` overestimates a digitized boundary's true length, so this caps around 0.90–0.98 even for a true disk and is size-dependent (0.977 at r=10px, 0.905 at r=160px, measured). Kept for backward compatibility with existing QC thresholds. Prefer `axon_circularity_crofton` for anything size-comparison-sensitive |
| `axon_circularity_crofton` | axon | `4π·area / perimeter_crofton²` | Shape roundness, unbiased estimator | ≈1.00 across the same size range where legacy `circularity` drifts from 0.98 to 0.90. Added alongside (not replacing) the legacy column |
| `axon_shape_irregularity` | axon | `1 / circularity` | Inverse of legacy circularity; `≥1`, increases with irregularity | Inherits legacy `circularity`'s size bias |
| `axon_shape_irregularity_crofton` | axon | `1 / axon_circularity_crofton` | Inverse of unbiased circularity | Preferred for cross-image/cross-size comparisons |
| `centroid_x_um`, `centroid_y_um` | axon | watershed centroid × `pixel_size_um` | Axon centroid position, calibrated | |
| `centroid_x_px`, `centroid_y_px` | axon | watershed centroid | Axon centroid position, raw pixels | |
| `nearest_neighbor_um` | axon | distance to nearest other axon centroid, this image | Local packing density proxy | `NaN` when the image has fewer than 2 axons |

---

## 3. Mitochondrial metrics — legacy per-axon aggregates (`axons.csv`)

Computed identically regardless of `--mito-assignment` (the underlying mitochondrion
list differs by mode, not the aggregation formulas). `ddof=0` for std/CV here — see
§4 for the Phase 3b metrics, which use `ddof=1`.

| metric_name | level | formula | interpretation | notes |
|---|---|---|---|---|
| `mito_count` | axon | count of mitochondrial objects assigned to this axon | | `0` is a real value, not missing data |
| `mito_area_um2` | axon | sum of assigned mitochondrial pixel areas × `pixel_size_um²` | Total mitochondrial area in this axon | Same quantity as `mito_total_area_um2` (§4), computed independently via the legacy code path — kept for backward compatibility |
| `mito_density_per_um2` | axon | `mito_count / axon_area_um2` | Mitochondria per unit axon area | `NaN` if `axon_area_um2 = 0` |
| `mito_mean_circularity` | axon | mean of `4π·area/perimeter²` over mitochondria | Mean roundness of mitochondria in this axon | `NaN` at `mito_count = 0` |
| `mito_mean_form_factor` | axon | mean of `perimeter²/(4π·area)` over mitochondria | Inverse of circularity; mitochondrial elongation/irregularity proxy | `NaN` at `mito_count = 0` |
| `mito_mean_feret_um` | axon | mean max Feret diameter over mitochondria | Mean mitochondrial "length" | `NaN` at `mito_count = 0` |
| `mito_std_area_um2` | axon | `std(mito areas, ddof=0)` | Size heterogeneity of mitochondria in this axon | `0.0` (not `NaN`) at `mito_count = 1`; `NaN` at `mito_count = 0` |
| `mito_max_area_um2` | axon | largest single mitochondrion's area | | `NaN` at `mito_count = 0` |
| `mito_cv_area` | axon | `mito_std_area_um2 / mean(areas)` | Coefficient of variation of mitochondrial size | `NaN` at `mito_count ≤ 1` |
| `mito_area_skewness` | axon | `scipy.stats.skew` of mitochondrial areas | Distribution asymmetry | `NaN` at `mito_count ≤ 2` (undefined below 3 points) |
| `mito_mean_dist_centroid_um` | axon | mean distance of mitochondrial centroids from the axon centroid | Central tendency of mitochondrial position | `NaN` at `mito_count = 0` |
| `mito_std_dist_centroid_um` | axon | `std` of those distances, `ddof=0` | Spread of mitochondrial position | `0.0` at `mito_count = 1`; `NaN` at `mito_count = 0` |

## 4. Mitochondrial burden & shape (`axons.csv`, Phase 3b)

Empty-mitochondria contract: count-based and ratio-based quantities are well-defined
zeros (`0 / positive = 0`); shape aggregates are `NaN` (undefined, not zero).

| metric_name | level | formula | interpretation | notes |
|---|---|---|---|---|
| `mito_total_area_um2` | axon | sum of mitochondrial areas | Total mitochondrial area | `0.0` at `mito_count = 0` (not `NaN`) |
| `mito_occupancy_ratio` | axon | `mito_total_area_um2 / axon_area_um2` | Fraction of axoplasm occupied by mitochondria | Mito-inclusive `axon_area_um2` under `mito_hole_handling=fill` (default) means this is a true areal fraction ∈ [0,1] |
| `mito_mean_area_um2` | axon | mean mitochondrion area | | `NaN` at `mito_count = 0` |
| `mito_median_area_um2` | axon | median mitochondrion area | Less sensitive to outlier mitochondria than the mean | `NaN` at `mito_count = 0` |
| `mito_area_iqr` | axon | `Q75 − Q25` of mitochondrial areas | Spread of mitochondrial size, robust to outliers | `NaN` at `mito_count = 0` |
| `mito_fragmentation_index` | axon | `mito_count / mito_total_area_um2` | Higher = same total mitochondrial area split into more, smaller pieces | `NaN` at `mito_count = 0` (never `inf`) |
| `normalized_mito_load` | axon | `mito_total_area_um2 / fiber_area_um2` | Mitochondrial burden relative to the whole fiber, not just axoplasm | Differs from `mito_occupancy_ratio`'s denominator |
| `mito_per_myelin` | axon | `mito_total_area_um2 / myelin_area_um2` | Mitochondrial burden relative to myelin | `NaN` whenever `myelin_area_um2` is `NaN` — **`NaN` for every pathological axon by design** (F4); do not treat as missing data to "fix" |
| `mito_mean_aspect_ratio` | axon | mean of `major_axis_length / minor_axis_length` | Mean mitochondrial elongation | `NaN` at `mito_count = 0`; guarded against `minor_axis_length = 0` |
| `mito_std_aspect_ratio` | axon | `std(aspect ratios, ddof=1)` | Spread of mitochondrial elongation | `NaN` at `mito_count ≤ 1` |
| `mito_mean_solidity` | axon | mean mitochondrial solidity | Mean mitochondrial boundary regularity | `NaN` at `mito_count = 0` |
| `mito_std_solidity` | axon | `std(solidity, ddof=1)` | | `NaN` at `mito_count ≤ 1` |
| `mito_mean_eccentricity` | axon | mean mitochondrial eccentricity | | `NaN` at `mito_count = 0` |

## 5. Mitochondrial spatial distribution (`axons.csv`, Phase 3b)

| metric_name | level | formula | interpretation | notes |
|---|---|---|---|---|
| `mito_peripheralization_index` | axon | `mito_mean_dist_centroid_um / equivalent_axon_radius_um` | ~0 = mitochondria central, ~1 = near the axon boundary | Values above 1 can occur for irregular axons (equivalent radius is approximate) — flag, don't discard |
| `mito_mean_nn_distance_um` | axon | mean nearest-neighbor distance among mitochondrial centroids | Observed mitochondrial spacing | `NaN` at `mito_count < 2` (undefined for 0–1 points) |
| `mito_clustering_index` | axon | `expected_nn_distance / observed_nn_distance` under a 2D Poisson (CSR) approximation | ~1 = random; >1 = clustered (observed NN smaller than expected); <1 = dispersed | `NaN` at `mito_count < 3` — the CSR approximation is explicitly unstable below that, even though `mito_mean_nn_distance_um` may already be available at count=2 |

---

## 6. Per-mitochondrion table (`mitochondria_metrics.csv`, `--write-mito-csv`)

One row per real mitochondrion — only produced under `--mito-assignment centroid`
or `overlap` (a per-fiber-crop fragment under `legacy` is not a real, whole
mitochondrion, so this table is empty in that mode). See `docs/known_issues.md` F8.

| metric_name | level | formula | interpretation | notes |
|---|---|---|---|---|
| `mito_id` | mitochondrion | global watershed label | Unique within one image | |
| `parent_axon_id` | mitochondrion | assigned axon's `axon_id` | Which axon this mitochondrion was attributed to | |
| `mito_area_um2` | mitochondrion | pixel area × `pixel_size_um²` | | |
| `mito_perimeter_um` | mitochondrion | perimeter × `pixel_size_um` | | |
| `mito_aspect_ratio` | mitochondrion | `major_axis_length / minor_axis_length` | Elongation | `NaN` if `minor_axis_length = 0` |
| `mito_solidity` | mitochondrion | `area / convex_hull_area` | Boundary regularity | |
| `mito_eccentricity` | mitochondrion | ellipse-fit eccentricity | | |
| `mito_circularity` | mitochondrion | `4π·area / perimeter²` | | Same digital-boundary bias as axon `circularity` (§2) |
| `mito_feret_um` | mitochondrion | max Feret diameter × `pixel_size_um` | Maximum caliper length | |
| `mito_centroid_x_px`, `mito_centroid_y_px` | mitochondrion | global pixel centroid | | |
| `mito_dist_to_axon_center_um` | mitochondrion | distance from this mitochondrion's centroid to its parent axon's centroid | | |

---

## 7. Image-level summary (`image_summary.csv`)

One row per image. `mean_*` columns are computed over the **QC-valid** subset by
default (identical to the `_all` counterpart when `--exclude-qc-failed` is off, since
nothing is excluded then); every `mean_*`/`image_*` column also has an `_all` twin
computed over every axon regardless of exclusion (see §7.3).

### 7.1 Legacy image aggregates

| metric_name | formula | interpretation | notes |
|---|---|---|---|
| `n_fibers` | count of detected axons | | |
| `fov_area_um2` | image height × width × `pixel_size_um²` | Field-of-view area | |
| `fiber_density_per_mm2` | `n_fibers / fov_area_um2 × 1e6` | Axons per mm² | |
| `mean_axon_area_um2` | mean of `axon_area_um2` | | |
| `mean_circularity` | mean of `circularity` | | Inherits legacy circularity's size bias — see §2 |
| `mean_convexity` | mean of `convexity` | | |
| `mean_avf` | mean of `axon_vol_fraction` | | |
| `mean_g_ratio` | mean of `g_ratio` | | `NaN` in pathological mode |
| `mean_myelin_thickness_um` | mean of `myelin_thickness_um` | | `NaN` in pathological mode |
| `mean_mvf` | mean of `myelin_vol_fraction` | | `NaN` in pathological mode |
| `mito_outside_area_um2` | mitochondrial pixel area outside axoplasm+myelin | Segmentation-noise / unassigned mitochondria indicator | Under `mito_hole_handling=legacy` this is structurally always the *entire* mitochondrial mask (see `docs/known_issues.md` F2) |
| `mito_outside_frac` | `mito_outside_area_um2 / total_mito_area` | Fraction of all mitochondrial pixels not attributed to any fiber | **Structurally always `1.0` under `--mito-hole-handling legacy`** — meaningful only under the default `fill` |
| `total_myelin_area_um2` | sum of myelin pixels in the whole mask | | |
| `myelin_area_fraction_of_fov` | `total_myelin_area_um2 / fov_area_um2` | | Feeds `image_demyelination_index` |
| `myelin_component_count` | connected myelin components in the mask | Fragmentation of the myelin sheath at the image level | |

### 7.2 Mitochondria-to-axon assignment diagnostics (Phase 3a)

| metric_name | formula | interpretation | notes |
|---|---|---|---|
| `n_mito_assigned` | count of mitochondria attributed to some axon | | `NaN` under `--mito-assignment legacy` (no global assignment computed in that mode) |
| `n_mito_unassigned` | total distinct mitochondria − `n_mito_assigned` | Mitochondria with no axon contact — check segmentation quality if this is large | `NaN` under `--mito-assignment legacy` |

### 7.3 Distribution statistics (Phase 4)

For each of `g_ratio`, `axon_area_um2`, `mito_density_per_um2`, `mito_occupancy_ratio`,
seven columns named `image_{stat}_{variable}`:

| stat | formula |
|---|---|
| `mean` | arithmetic mean |
| `std` | standard deviation, `ddof=1` |
| `cv` | `std / mean` |
| `median` | median |
| `min` | minimum |
| `max` | maximum |
| `iqr` | `Q75 − Q25` |

e.g. `image_cv_g_ratio`, `image_median_axon_area_um2`, `image_iqr_mito_occupancy_ratio`.
All `NaN` on empty input; never raises. A single valid axon yields `std`/`cv` = `NaN`
(`ddof=1` is undefined for n=1) but `iqr = 0`.

### 7.4 MVF and demyelination index (Phase 4)

| metric_name | formula | interpretation | notes |
|---|---|---|---|
| `image_mvf` | `Σ myelin_area_um2 / Σ fiber_area_um2` over valid axons | Image-level myelin volume fraction | `NaN` if every axon's `myelin_area_um2` is `NaN` (e.g. an all-pathological image) |
| `image_demyelination_index` | `1 − (myelin_area_fraction_of_fov / reference)`, clipped to [0,1] | Demyelination relative to a reference (normal-tissue) myelin fraction | **`NaN` unless `--demyelination-reference` is explicitly supplied.** No default — this is a biological calibration choice the tool cannot derive on its own (see `docs/known_issues.md` F4) |

### 7.5 QC summary (Phase 5) — see `docs/qc.md` for the flags themselves

| metric_name | formula | interpretation | notes |
|---|---|---|---|
| `image_n_axons_total` | `n_fibers` | | |
| `image_n_axons_valid` | axons with `excluded_from_analysis = False` | | Equals total when `--exclude-qc-failed` is off |
| `image_n_axons_excluded` | `total − valid` | | Always `0` unless `--exclude-qc-failed` is set |
| `image_percent_flagged_axons` | `100 × (fraction of axons with any qc_* flag True)` | | Flags ≠ exclusions — this counts flags regardless of exclusion setting |
| `qc_high_exclusion_rate` | `True` if `excluded / total > 0.3` | Per-image warning: over 30% of axons excluded | |
| `image_noncanonical_mask_frac` | fraction of mask pixels outside `{0, myelin_val, mito_val, axoplasm_val}` | Mask quality indicator | See `docs/known_issues.md` F6 |
| `qc_mask_noncanonical` | `True` if `image_noncanonical_mask_frac > 0.001` | | |

### 7.6 `_all`-suffixed unfiltered aggregates

Every column in §7.1 (the `mean_*` block) and §7.3 (`image_*` distribution stats) has
an `_all` counterpart (`mean_axon_area_um2_all`, `image_cv_g_ratio_all`, ...) computed
over **every** axon regardless of `excluded_from_analysis`. Always present, identical
to the non-`_all` version when `--exclude-qc-failed` is off, and the only way to see
the unfiltered picture when it's on.

---

## 8. `group_metrics.csv` (`--write-group-csv`)

One row per `group` value present in `image_summary.csv`: `n_images`, plus `{metric}_mean`
and `{metric}_std` (`ddof=1`) for a fixed set of key image-level metrics (mean g-ratio,
mean axon area, mean AVF/MVF, mito occupancy/density distribution stats, `image_mvf`,
`image_demyelination_index`, `mito_outside_frac`). See `temptation/summaries.py:summarize_groups`
for the exact list.

---

## 9. QC flags and exclusion columns

`qc_*` boolean flags, `excluded_from_analysis`, and `exclusion_reason` are documented
separately in **`docs/qc.md`** — they follow their own rules (which flags can exclude,
which never can, and why) that don't fit this table's format cleanly.
