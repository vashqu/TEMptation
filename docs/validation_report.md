# Validation report — normal vs. pathological, this repository's dataset

**Sanity, not significance.** Everything below is descriptive: ranges, means, and
whether two specific pre-registered checks hold. No hypothesis test, no p-value, no
ROC curve, no classifier accuracy figure appears anywhere in this document, on
purpose — see CLAUDE.md §9/§20 and `docs/pathology_score.md`. The corresponding
automated checks live in `tests/test_biological_sanity.py`, all passing as of this
report.

**Run**: full `normal_data`/`pathological_data` trees (`--input-root .. --groups
normal_data:normal pathological_data:pathological --pixel-um 0.00524`), default
`SegmentationConfig`/`QCThresholds` (i.e. `--mito-hole-handling fill`,
`--mito-assignment centroid`). 100 image pairs, 0 skipped, 0 processing errors:
94 normal images / 1140 axons, 6 pathological images / 25 axons.

---

## 1. CLAUDE.md §13 invariants — all hold on the normal group

| invariant | result |
|---|---|
| `0 ≤ g_ratio ≤ 1` | holds for all 1140/1140 axons (range 0.287–0.875) |
| `fiber_area_um2 > axon_area_um2` | holds for all 1140 axons |
| `myelin_area_um2 = fiber_area_um2 − axon_area_um2` | max residual `8.9e-16` (floating-point identity) |
| `0 ≤ mito_occupancy_ratio ≤ 1` | holds (range 0.000–0.853) |
| `0 ≤ normalized_mito_load ≤ 1` | holds (range 0.000–0.488) |
| `axon_shape_irregularity ≥ 1` | holds (min 1.104); Crofton variant also holds (min 1.016) |
| `circularity ≤ 1` | holds (max 0.906 — see `docs/known_issues.md` F5 for why this doesn't reach 1.0); Crofton variant holds under a 1.05 tolerance (max 0.984) |

No violations anywhere. `tests/test_biological_sanity.py` pins all seven as
regression tests.

## 2. Pre-registered baseline: `mean_g_ratio`

CLAUDE.md's original audit estimated `mean_g_ratio ≈ 0.63` from a 4-image subset
(`normal_data/162-165`). Confirmed at full scale (n=1140 axons across 94 images):

```
mean_g_ratio = 0.679
```

Comfortably inside the documented acceptable range `[0.55, 0.75]`. (The 4-image
subset alone gives 0.655 — closer to the original estimate, and also inside range;
the difference is just sample composition, not a discrepancy worth chasing.)

## 3. Documented degeneracies on the pathological group — confirmed present

Per `docs/known_issues.md` F4, myelin is fully detached in this group, and three
things should be true as a direct consequence. All three confirmed on all 25
pathological axons:

- `axon_vol_fraction` — exactly `{1.0}`, no other value present.
- `myelin_area_um2`, `g_ratio`, `myelin_vol_fraction` — entirely `NaN`.
- `qc_no_myelin` — `True` for every axon (the permanently-non-excludable flag that
  keeps this group from being silently deleted under `--exclude-qc-failed`; see
  `docs/qc.md`).

These are pinned as regression tests too — a future segmentation change that
"accidentally" produces a non-degenerate value here should fail loudly, not pass
silently, since it would mean this dataset's myelin is no longer being treated as
detached, which would itself need `docs/known_issues.md` updated deliberately, not
discovered by surprise.

## 4. `image_demyelination_index` — strictly separates the groups on this dataset

`image_demyelination_index` needs an externally-supplied reference (CLAUDE.md §11
A4 deliberately gives `--demyelination-reference` no default for real analyses —
this is a biological calibration choice, not something the tool should guess).
For this validation run only, the reference was the **normal group's own median**
`myelin_area_fraction_of_fov` (0.264):

```
normal:        mean 0.117   range [0.000, 0.727]
pathological:  mean 0.950   range [0.884, 1.000]
```

Pathological's minimum (0.884) exceeds normal's maximum (0.727) — strict
separation on this dataset, matching the blueprint's pre-registered check.
**This reference value is a convenience choice for this validation run, not a
general-purpose default** — a different imaging setup or staining protocol would
need its own reference derived the same way (median `myelin_area_fraction_of_fov`
over a known-normal reference set), which is exactly why the CLI flag has no
built-in default.

## 5. Pathology-score component distributions — why three stay blocked

`docs/pathology_score.md` lists three score components with no calibrated
threshold (`low_myelin_fraction`, `high_mito_fragmentation`, `high_mito_occupancy`).
Real distributions, inspected here specifically to check whether this validation
run supplies enough basis to unblock any of them:

**`mito_fragmentation_index`** (axon-level) — heavily right-skewed, and only
defined where `mito_count > 0` (592/1140 normal axons; 10/25 pathological axons):

```
normal:        median 23.2   p90 267.6   max 1655.4   (n=592 valid)
pathological:  median 15.6                            (n=10 valid)
```

A three-order-of-magnitude range with the bulk of mass far from the tail makes any
single cutoff highly sensitive to exactly where you draw it, and the small
valid-n subset (especially pathological's 10/25) leaves little basis for choosing
one. **Stays blocked.**

**`mito_occupancy_ratio`** (axon-level) — the more striking finding:

```
normal:        mean 0.069   median 0.008   p90 0.222   max 0.853
pathological:  mean 0.042   median 0.000   max 0.230   (n=25)
```

Pathological's mean is **lower**, not higher, than normal's — the opposite
direction a "high occupancy indicates pathology" component would assume. Setting a
threshold here based on an assumption this dataset doesn't support would encode a
wrong prior into the score. **Stays blocked** — and this finding should inform
*whether* this component belongs in the score at all, not just where to draw its
line, whenever a biologist next reviews `docs/pathology_score.md`.

**`low_myelin_fraction`** — the axon-level `myelin_vol_fraction` is structurally
`NaN` for every pathological axon (same F4 degeneracy as §3), so it can only ever
flag *normal* axons with unusually low individual myelin fraction — not useful for
a cross-group indicator without the same externally-supplied-reference treatment
`image_demyelination_index` already needed (§4). **Stays blocked** for the same
underlying reason, not a new one.

## 6. Recommendation

No pathology-score component gained a defensible threshold from this pass. The
three already-usable components (`g_ratio_high`, `g_ratio_low`, `low_circularity`,
which reuse `QCThresholds`' existing calibrated values) are the only ones with real
thresholds behind them — see `docs/pathology_score.md` for whether/when a v0 score
built from just those three is worth doing; that's a product decision (does a
score built from 3-of-6 intended components, none of them mitochondrial, say
anything useful?), not a data question this report can settle by itself.
