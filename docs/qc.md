# Quality control: flags and exclusion

Two separate concerns, computed by `temptation/qc.py`:

- **Flagged** — every `qc_*` boolean column is always computed for every axon, in
  every run. A flag is an observation, not an action. Flags never remove data.
- **Excluded** — only happens when you explicitly pass `--exclude-qc-failed` (CLI)
  or tick **"Exclude QC-failed objects from image summaries"** (GUI, Setup → QC).
  Excluded axons are *marked*, never dropped — `axons.csv` has the same row count
  whether or not `--exclude-qc-failed` is set. This is what makes it possible to
  compute both a raw and a QC-filtered image summary from the same table (the
  `_all`-suffixed columns in `image_summary.csv` — see `docs/metrics.md` §7.6).

## The flags

| flag | condition | can exclude? |
|---|---|---|
| `qc_g_ratio_low` | `g_ratio < g_ratio_min` (default 0.3) | yes |
| `qc_g_ratio_high` | `g_ratio > g_ratio_max` (default 0.95) | yes |
| `qc_extreme_g_ratio` | `qc_g_ratio_low OR qc_g_ratio_high` | yes |
| `qc_low_circularity` | `circularity < axon_circularity_min` (default 0.4) | yes |
| `qc_invalid_area_relation` | `axon_area_um2 > fiber_area_um2` (strict `>`) | yes |
| `qc_low_myelin_area` | `myelin_area_um2 <= min_myelin_area_um2` | yes |
| `qc_low_axon_area` | `axon_area_um2 <= min_axon_area_um2` | yes |
| `qc_low_fiber_area` | `fiber_area_um2 <= min_fiber_area_um2` | yes |
| `qc_no_mitochondria` | `mito_count == 0` | **never** |
| `qc_low_mito_count_clustering` | `mito_count < min_mito_count_for_clustering` (default 3) | **never** |
| `qc_no_myelin` | `myelin_area_um2` is `NaN` | **never** |

Thresholds live in `temptation.config.QCThresholds` and are set via CLI flags
(`--g-ratio-min`, `--g-ratio-max`, `--axon-circularity-min`, `--min-myelin-area-um2`,
`--min-axon-area-um2`, `--min-fiber-area-um2`, `--min-mito-count-for-clustering`) or the
GUI's Setup → QC tab, which reflects every `QCThresholds` field automatically.

**A threshold of `None` (the default for the three `min_*_area_um2` fields) means "not
calibrated for this dataset" — the flag stays all-`False`, not all-`True`.** These three
need a value specific to your imaging setup before they mean anything; don't set them to
`0` expecting that to be a no-op, it silently behaves the same as `None` for the `<=`
comparison used here, but `None` is the explicit, documented "not yet calibrated" state.

## The trap: three flags that must never exclude

`qc_no_myelin`, `qc_no_mitochondria`, and `qc_low_mito_count_clustering` are permanently
excluded from `EXCLUDABLE_FLAGS`, regardless of `--exclude-qc-failed`. This is deliberate,
not an oversight:

- `qc_no_myelin` fires on **every single pathological axon** — myelin is detached, not a
  data-quality problem (see `docs/known_issues.md` F4). Making this excludable would
  silently delete the entire pathological group under `--exclude-qc-failed`.
- `qc_low_mito_count_clustering` fires whenever fewer than 3 mitochondria are present,
  which is the *common* case for a typical axon crop, not rare or defective.
- `qc_no_mitochondria` (0 mitochondria) is frequently a real biological observation,
  not a segmentation failure.

If you need to filter these cases out for a specific analysis, do it downstream from the
CSV — don't ask this tool to exclude them, since that would conflate "this axon has an
uninteresting/undesirable value" with "this axon's segmentation is untrustworthy," which
is what exclusion is for.

## `axon_area_um2 > fiber_area_um2`: strict `>`, not `>=`

CLAUDE.md's original spec used `>=`. In practice `axon_area_um2 == fiber_area_um2` is a
**structural invariant**, not an anomaly, under the default `--assign-detached-myelin
none`: it's exactly what happens whenever zero myelin was captured for a fiber, which
includes every pathological axon and any "normal"-mode image whose myelin has detached
(F3/F4). That case is already tracked (and deliberately non-excludable) via `qc_no_myelin`.
Using `>=` here made `qc_invalid_area_relation` redundant with `qc_no_myelin` in exactly
the one case that must never exclude — except this flag *was* excludable, so it silently
excluded 100% of a pathological image's axons under `--exclude-qc-failed` before this was
caught. Strict `>` keeps the flag meaningful as a genuine geometric-impossibility check,
which can actually occur under `--assign-detached-myelin nearest` (where `fiber_area_px`
is computed independently rather than as a superset of the axon).

## Exclusion reasons

When `--exclude-qc-failed` is set, `exclusion_reason` is a `;`-joined, alphabetically
sorted list of short tokens for every *excludable* flag that fired:

```
g_ratio_high;low_circularity
```

`excluded_from_analysis = True` exactly when `exclusion_reason` is non-empty. Both
columns are always present in `axons.csv`, regardless of `--exclude-qc-failed` —
`excluded_from_analysis` is just `False`/`""` for every row when the flag is off.

`run_manifest.json`'s `exclusion_counts` field is a per-token tally across the whole
batch (an axon excluded for two reasons counts toward both tokens), so an
exclusion-heavy run is auditable at a glance without scanning every row.

## `qc_report.csv` (`--write-qc-report`)

A slim view of `axons.csv`: identifiers (`image_id`, `group`, `axon_id`) + every
`qc_*` flag + `excluded_from_analysis`/`exclusion_reason`, with none of the metric
columns. Useful for a quick pass/fail scan without opening the full table.

## `image_percent_flagged_axons` vs `qc_high_exclusion_rate`

These answer different questions and are easy to conflate:

- `image_percent_flagged_axons` — what fraction of this image's axons have **any**
  `qc_*` flag `True` (including the three never-excludable ones). Computed regardless
  of `--exclude-qc-failed`.
- `qc_high_exclusion_rate` — `True` if more than 30% of this image's axons were
  **excluded**. Always `False` unless `--exclude-qc-failed` is set (nothing is ever
  excluded otherwise).

An image can have a very high flagged percentage (e.g. an all-pathological image,
where `qc_no_myelin` fires on every axon) while `qc_high_exclusion_rate` stays `False`
the entire time, because that flag can never exclude. That's expected, not a bug —
see the trap above.

## GUI live preview

The Setup → QC tab's live preview (blueprint-era feature) recomputes flag counts
against the last completed run's raw per-axon metrics whenever you edit a threshold —
no re-segmentation, and it survives later parameter edits without needing a fresh Run.
It calls `temptation.qc.compute_qc_flags` directly, so the preview can never drift from
what an actual Run would produce.
