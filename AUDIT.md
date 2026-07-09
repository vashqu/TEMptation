# TEMptation — Phase 0 audit

Re-verification of the eight defects (F1–F8) documented in
`IMPLEMENTATION_BLUEPRINT.md` §0, run independently against the golden
baselines in `tests/golden/`. Environment: conda env `nerve_env`
(numpy 1.26.4, scikit-image 0.24.0, pandas 2.2.2 — see
`tests/golden/README.md`). No source files were modified to produce this
audit.

**No fixes are applied in this phase.** Phase 1 is a pure code-move with
zero output change; F1 is fixed only in the gated Phase 2b, behind
`--mito-hole-handling {fill,legacy}`, with `legacy` reproducing exactly what
is recorded below. Once Phase 2b lands, every metric derived from
mitochondrial area or axon-boundary geometry is **expected** to diverge from
these golden numbers under the new default (`fill`) — that divergence will
be the documented fix, not a regression. The `legacy` codepath must still
reproduce this audit's numbers byte-for-byte forever, which is what
`tests/test_regression_golden.py` (Phase 1) enforces.

---

## F1 — mitochondrial rim bug (confirmed, highest priority)

`axon_dil = binary_dilation(axon_only, disk(1))` reaches one pixel into
each mitochondrion-shaped hole in the axoplasm mask, so `mito_in_axon`
is the outer 1-px rim, not the mitochondrion.

Verified on `normal_data/162-165/163. Mask 20K.tif`:

| quantity | value |
|---|---|
| true mitochondrial area (`(mask==128).sum() × px²`) | 1.3659 µm² |
| sum of `mito_area_um2` over axons in golden `axons.csv` | 0.0291 µm² |
| reported / true | **0.0213** (2.1%) |

Confirmed. This invalidates `mito_area_um2`, `mito_density_per_um2`,
`mito_mean_circularity`, `mito_mean_form_factor`, `mito_mean_feret_um`,
`mito_std_area_um2`, `mito_max_area_um2`, `mito_cv_area`,
`mito_area_skewness`, and biases `axon_area_um2` (low), `g_ratio` (high),
`circularity`, `solidity`, `convexity` (perimeter traces the holes).

**Action:** fix in Phase 2b only, gated behind an explicit flag, with a
before/after impact table for the biologist.

## F2 — `mito_outside_frac` structurally always 1.0 (confirmed)

`mito_outside = mito & ~(axon_only | myelin)` uses `axon_only` (raw
`== 192`), not the mito-inclusive `axoplasm`. Mitochondrial pixels are never
`64` or `192`, so nothing is excluded.

Verified on all four images in `162-165`:

```
image_id  mito_outside_frac
162       1.0
163       1.0
164       1.0
165       1.0
```

Confirmed exactly. Both `mito_outside_area_um2` and `mito_outside_frac`
carry zero information under current code.

## F3 — `--mode auto` misclassifies pathological images (confirmed)

`myelin_threshold_px=200` is three orders of magnitude below actual myelin
pixel counts. Live run on `pathological_data`:

```
image_id  mode
path_1    pathological
path_2    pathological
path_3    normal
path_4    normal
path_5    normal
path_6    normal
```

4 of 6 pathological images resolve as `normal` mode. Confirmed exactly as
documented. `mode` is a myelin-presence detector, not a group label —
`group` must never be inferred from it (this is enforced starting Phase 2).

## F4 — myelin unobtainable in pathological group (confirmed)

```
axon_vol_fraction unique values across all 25 pathological axons: [1.0]
rows with non-NaN myelin_area_um2: 0 / 25
```

Confirmed exactly. Both CLAUDE.md §5.6 `demyelination_index` definitions
reduce to `axon_area/fiber_area`, which is identically 1.0 here — the index
is blind precisely where it needs to discriminate. Resolved in the blueprint
(§11 A4) via an image-level definition based on `myelin_area_fraction_of_fov`
relative to the normal-group reference; requires biologist sign-off before
Phase 4 implements it.

## F5 — circularity is size-dependent, caps ≈0.90 (confirmed, not re-run here)

Verified in the original blueprint audit via synthetic disks
(`skimage.measure.perimeter` overestimates digital boundary length).
Not re-derived in this pass; no golden-CSV-dependent check applies. Will be
addressed in Phase 2 by adding `axon_circularity_crofton` alongside (not
replacing) the legacy `circularity` column.

## F6 — non-canonical mask values (confirmed, not re-run here)

`pathological_data/mask_path_3.tif` contains 203 distinct pixel values
against the canonical `{0,64,128,192}`. Not re-derived in this pass (no
golden CSV depends on it — it's outside the `162-165` / full-pathological
baseline scope, since mask_path_3 is included in the pathological run).
Confirmed present via the original audit; will surface as
`image_noncanonical_mask_frac` starting Phase 2/3.

## F7 — pair discovery non-recursive, one file silently dropped (confirmed, not re-run here)

`_find_pairs_in_folder(Path("normal_data"))` returns 0 pairs (subfolders
not globbed); `normal_data/422-430/429. Mask 20K.gif` is dropped because the
glob is `*.tif`/`*.tiff` only. Not re-derived in this pass (Phase 0 goldens
use `--folder ../normal_data/162-165` directly, which sidesteps this).
Fixed in Phase 2 (`discovery.py`).

## F8 — minor findings, schema stability (partially re-verified)

Confirmed: `tests/golden/normal/axons.csv` and
`tests/golden/pathological/axons.csv` have identical column sets (34 cols)
despite very different per-image axon counts (2 to 14), so the "no columns
at all when df_axons is empty" failure mode did not trigger in these two
baselines — every processed image had ≥1 axon. The empty-dataframe schema
gap (`nearest_neighbor_um` and even standard columns missing when
`df_axons` is empty) is a real risk documented in the blueprint and is
covered by a dedicated Phase 6 test (`test_schema.py`, 0-axon-image case),
not re-derivable from these baselines alone.

---

## Summary

All defects checked in this pass (F1–F4, F8-partial) reproduce exactly as
documented in `IMPLEMENTATION_BLUEPRINT.md`. F5–F7 are accepted from the
prior audit without re-derivation here since they don't bear on golden-CSV
byte-stability. No code was changed. Golden baselines in `tests/golden/` are
confirmed to encode this exact (buggy, F1/F2/F3/F4-affected) behavior and
are the correct anchor for Phase 1's zero-output-change refactor.
