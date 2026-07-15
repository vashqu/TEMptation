# Morphometric pathology score (planned — not yet implemented)

> **Status: not built.** No `morphometric_pathology_score` column exists in any CSV
> today, no `--pathology-score` CLI flag exists, and no `pathology_score.py` module
> exists. `temptation.config.PathologyScoreConfig`/`ComponentSpec` are scaffolding
> dataclasses only — they aren't read by anything yet. This document describes the
> plan so a future implementer (human or model) has the design decided in advance,
> and so this feature is never accidentally presented as more validated than it is.

## The one thing that must never change about this feature

> **This is not a diagnostic score. It is a morphometric abnormality score based on
> predefined mask-derived features.**

That sentence must appear verbatim in every surface this feature ever reaches:
the CSV column's documentation (this file), the GUI (wherever the score is displayed
or enabled), the CLI help text for whatever flag enables it, and `run_manifest.json`'s
record of the run. If you implement this feature and any of those four surfaces is
missing the sentence, the implementation is incomplete regardless of whether the
math is correct.

## Why this isn't built yet

CLAUDE.md's own priority order (§11) places this last, after QC, image-level
summaries, and spatial metrics — deliberately: a pathology score is only as
trustworthy as the thresholds and weights that go into it, and several of those
inputs need dataset-specific calibration this tool cannot supply on its own (see
`docs/qc.md`'s discussion of `None`-valued QC thresholds, which is the same
underlying issue). Phase 9 (biological validation, not yet started as of this
writing) is where real-data distributions get inspected before any threshold here
is set to a non-`None` value.

## Planned design (rule-based, not a classifier)

**Do not train a classifier.** The score must stay rule-based and interpretable —
a weighted sum of independently-inspectable abnormality indicators, each of which
a biologist can look at and agree or disagree with individually. If a future
implementer is tempted to fit weights so that `normal` and `pathological` groups
separate better, that is hand-fitting a classifier to two labels and is explicitly
out of scope; see `docs/known_issues.md`'s note on Phase 9's validation risk.

Candidate components (from `temptation.config.ComponentSpec`: a metric name, a
comparison operator, a threshold, and a weight):

| component | candidate threshold | candidate weight | status |
|---|---|---|---|
| high g-ratio | `g_ratio > 0.95` | 20 | threshold has a real default (`QCThresholds.g_ratio_max`) — usable |
| low g-ratio | `g_ratio < 0.3` | 10 | threshold has a real default (`QCThresholds.g_ratio_min`) — usable |
| low myelin fraction | *(none set)* | 20 | **blocked** — needs a calibrated reference, same problem as `image_demyelination_index` (`docs/metrics.md` §7.4) |
| high mitochondrial fragmentation | *(none set)* | 20 | **blocked** — needs empirical distribution inspection |
| high mitochondrial occupancy | *(none set)* | 15 | **blocked** — needs empirical distribution inspection |
| low axon circularity | `axon_circularity < 0.4` | 15 | threshold has a real default (`QCThresholds.axon_circularity_min`) — usable |

Any component whose threshold is `None` must **not** contribute to the score until
a biologist supplies a real value or an empirical distribution is inspected on real
data — silently defaulting an unset threshold to some guessed number would make the
score's output depend on an implementation detail nobody reviewed.

## Planned output columns

```
morphometric_pathology_score        (0-100)
morphometric_pathology_score_version   (e.g. "v0")
pathology_score_components          (which components fired, and their individual contribution)
```

`pathology_score_components` matters as much as the score itself — a `72` with no
way to see *why* is far less useful to a biologist than a `72` that breaks down into
"high g-ratio (+20), low circularity (+15), ...". This mirrors `exclusion_reason`'s
own design (`docs/qc.md`): always show the reasons, never just the verdict.

## Planned CLI/GUI surface

- CLI: an opt-in flag (e.g. `--pathology-score`), off by default. `PathologyScoreConfig.enabled`
  already exists in `temptation/config.py` for this.
- GUI: nothing today. When built, it should live in Setup (a way to enable it and see
  which components are active/blocked) and the Results dashboard (the score's
  distribution across groups — **descriptively**, not as a hypothesis test or
  classifier accuracy figure; see the note on Phase 9 below).

## What "done" looks like (Phase 9 gate)

Per CLAUDE.md §9/§20: implement the score only after QC thresholds and image-level
summaries have real data behind them, run it on both `normal_data` and
`pathological_data`, and **report the distribution descriptively — no hypothesis
test, no ROC curve, no accuracy number.** The temptation to tune weights until the
two groups separate is exactly the thing this whole document exists to prevent.
