# Golden baseline outputs

These CSVs are the frozen, pre-refactor outputs of `measure_nerve.py`. Every
later phase of `IMPLEMENTATION_BLUEPRINT.md` must reproduce them exactly
(within `rtol=1e-12`) whenever run with `--mito-hole-handling legacy` /
`--schema 1.0`, since Phase 1 introduces zero behavior change and later
phases must remain switchable back to legacy behavior.

## Environment (conda env `nerve_env`)

```
python        3.11
numpy         1.26.4
scipy         1.13.1
scikit-image  0.24.0
pandas        2.2.2
tifffile      2024.7.2
matplotlib    3.9.1
```

Interpreter used: `/opt/anaconda3/envs/nerve_env/bin/python`

## Commands used to generate these files

Run from `TEMptation/`:

```bash
python measure_nerve.py --folder ../normal_data/162-165 --pixel-um 0.00524 \
    --output-dir tests/golden/normal

python measure_nerve.py --folder ../pathological_data --pixel-um 0.00524 \
    --output-dir tests/golden/pathological
```

## Expected shapes

| file | rows | columns |
|---|---|---|
| `normal/axons.csv` | 47 | 34 |
| `normal/image_summary.csv` | 4 | 17 |
| `pathological/axons.csv` | 25 | 34 |
| `pathological/image_summary.csv` | 6 | 17 |

Per-image axon counts (normal `162-165`): 162→6, 163→14, 164→14, 165→13.
Per-image axon counts (pathological): path_1→7, path_2→2, path_3→5,
path_4→3, path_5→4, path_6→4.

`--mode auto` resolves `path_1`/`path_2` as `PATHOLOGICAL` and
`path_3`..`path_6` as `NORMAL` — this is expected (see `AUDIT.md` finding
F3) and is not a bug to fix in these baselines.

## Determinism

Verified: running the normal-group command twice into different output
directories produces byte-identical `axons.csv` and `image_summary.csv`.
The pipeline has no randomness (watershed is deterministic given fixed
inputs), so any future divergence in a regression test indicates a real
behavior change, not run-to-run noise.
