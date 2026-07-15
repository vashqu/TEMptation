"""Biological sanity checks (CLAUDE.md Sec 13, blueprint Phase 9): run
both real groups end to end and assert invariants that must hold if the
pipeline is behaving correctly, plus the documented *degeneracies* on
the pathological group (CLAUDE.md Sec 0 F4) so a future change that
silently "fixes" them fails loudly instead of quietly changing what
this tool measures.

Sanity, not significance -- no hypothesis test, no accuracy figure, no
group-separation claim beyond the two explicit, pre-registered checks
below (mean_g_ratio's documented range, and demyelination_index's
strict separation). See docs/validation_report.md for the full
descriptive writeup this test's numbers feed into.
"""

from pathlib import Path

import pytest

from temptation import discovery, pipeline
from temptation.config import QCThresholds, SegmentationConfig

DATA_ROOT = Path(__file__).parent.parent.parent
GROUP_MAP = {"normal_data": "normal", "pathological_data": "pathological"}


@pytest.fixture(scope="module")
def full_dataset_result():
    """Runs once for the whole module (not once per test) -- a full
    100-image batch takes on the order of a minute, and every test here
    reads from the same run."""
    pairs, skipped = discovery.find_groups(DATA_ROOT, GROUP_MAP)
    assert skipped == [], f"expected 0 skipped files, got {skipped}"
    result = pipeline.analyze_dataset(
        pairs, pixel_length_um=0.00524,
        seg_cfg=SegmentationConfig(), qc_thresholds=QCThresholds(),
    )
    assert result.errors == [], f"expected 0 processing errors, got {result.errors}"
    return result


@pytest.fixture(scope="module")
def normal_axons(full_dataset_result):
    df = full_dataset_result.df_axons
    return df[df["group"] == "normal"]


@pytest.fixture(scope="module")
def pathological_axons(full_dataset_result):
    df = full_dataset_result.df_axons
    return df[df["group"] == "pathological"]


# ---------------------------------------------------------------- CLAUDE.md Sec 13

def test_normal_g_ratio_in_unit_interval(normal_axons):
    g = normal_axons["g_ratio"].dropna()
    assert len(g) == len(normal_axons), "g_ratio should be fully populated for normal-mode axons"
    assert (g >= 0).all() and (g <= 1).all()


def test_normal_fiber_area_exceeds_axon_area(normal_axons):
    assert (normal_axons["fiber_area_um2"] > normal_axons["axon_area_um2"]).all()


def test_normal_myelin_equals_fiber_minus_axon(normal_axons):
    residual = (
        normal_axons["myelin_area_um2"]
        - (normal_axons["fiber_area_um2"] - normal_axons["axon_area_um2"])
    ).abs()
    assert residual.max() < 1e-9, f"max residual {residual.max()!r} exceeds float tolerance"


def test_normal_mito_occupancy_in_unit_interval(normal_axons):
    occ = normal_axons["mito_occupancy_ratio"].dropna()
    assert (occ >= 0).all() and (occ <= 1).all()


def test_normal_normalized_mito_load_in_unit_interval(normal_axons):
    load = normal_axons["normalized_mito_load"].dropna()
    assert (load >= 0).all() and (load <= 1).all()


def test_normal_axon_shape_irregularity_at_least_one(normal_axons):
    irregularity = normal_axons["axon_shape_irregularity"].dropna()
    assert (irregularity >= 1).all()
    irregularity_crofton = normal_axons["axon_shape_irregularity_crofton"].dropna()
    assert (irregularity_crofton >= 1).all()


def test_normal_circularity_at_most_one(normal_axons):
    """Legacy `circularity` is a biased estimator that never exceeds 1 in
    practice (it *undershoots*, per docs/known_issues.md F5) -- this
    assertion is the invariant CLAUDE.md Sec 13 asks for, not a claim
    that 1.0 means "perfect circle" (it doesn't; see docs/metrics.md)."""
    circularity = normal_axons["circularity"].dropna()
    assert (circularity <= 1).all()
    # axon_circularity_crofton is the unbiased estimator and can drift
    # slightly above 1.0 on a real (non-perfectly-circular, digitized)
    # boundary -- a documented tolerance, not a bug. 1.05 matches the
    # tolerance blueprint used for the same reason.
    circularity_crofton = normal_axons["axon_circularity_crofton"].dropna()
    assert (circularity_crofton <= 1.05).all()


def test_normal_mean_g_ratio_within_documented_range(normal_axons):
    """Pre-registered baseline check, not tuned after seeing the number:
    blueprint's original audit estimated ~0.63 from a 4-image subset.
    Confirmed here at n=1140 (the full normal_data group): 0.679, still
    comfortably inside [0.55, 0.75]. If a future segmentation change
    moves this outside that range, that's worth a human look before
    merging, not necessarily a bug."""
    mean_g_ratio = normal_axons["g_ratio"].mean()
    assert 0.55 <= mean_g_ratio <= 0.75, f"mean_g_ratio={mean_g_ratio!r} outside documented range"


# ---------------------------------------------------------------- documented degeneracies (F4)

def test_pathological_axon_vol_fraction_is_degenerate(pathological_axons):
    """CLAUDE.md Sec 0 F4: myelin is fully detached in this group, so
    axon_vol_fraction (axon_area/fiber_area) collapses to exactly 1.0
    for every axon -- fiber_area_px is defined to equal axon_area_px
    when watershed has no myelin channel to expand into. This is a
    documented, expected degeneracy, not something a future change
    should silently "fix" without updating docs/known_issues.md F4 and
    this test together."""
    assert sorted(pathological_axons["axon_vol_fraction"].unique()) == [1.0]


def test_pathological_myelin_is_entirely_nan(pathological_axons):
    assert pathological_axons["myelin_area_um2"].isna().all()
    assert pathological_axons["g_ratio"].isna().all()
    assert pathological_axons["myelin_vol_fraction"].isna().all()


def test_pathological_qc_no_myelin_fires_on_every_axon(pathological_axons):
    """The flip side of the degeneracy above: qc_no_myelin (docs/qc.md)
    must catch every one of these axons, since it's the mechanism that
    keeps them from being silently excluded under --exclude-qc-failed
    (qc_no_myelin is permanently non-excludable)."""
    assert pathological_axons["qc_no_myelin"].all()


# ---------------------------------------------------------------- demyelination_index

def test_demyelination_index_separates_groups_with_a_real_reference():
    """Not a hypothesis test -- a single descriptive check that the
    index, calibrated against a reference derived from the normal
    group's own data, produces non-overlapping ranges on this specific
    dataset. See docs/validation_report.md for the actual numbers and
    the caveat that this reference (median of the normal group's own
    myelin_area_fraction_of_fov) is a convenience choice for this
    validation run, not a general-purpose default -- CLAUDE.md Sec 11 A4
    deliberately left --demyelination-reference with no default for
    real analyses."""
    pairs, _ = discovery.find_groups(DATA_ROOT, GROUP_MAP)
    normal_pairs = [p for p in pairs if p["group"] == "normal"]
    prelim = pipeline.analyze_dataset(
        normal_pairs, pixel_length_um=0.00524,
        seg_cfg=SegmentationConfig(), qc_thresholds=QCThresholds(),
    )
    reference = prelim.df_image["myelin_area_fraction_of_fov"].median()

    result = pipeline.analyze_dataset(
        pairs, pixel_length_um=0.00524,
        seg_cfg=SegmentationConfig(), qc_thresholds=QCThresholds(),
        reference_myelin_fraction=reference,
    )
    df_image = result.df_image
    normal_di = df_image.loc[df_image["group"] == "normal", "image_demyelination_index"]
    patho_di = df_image.loc[df_image["group"] == "pathological", "image_demyelination_index"]

    assert normal_di.notna().all() and patho_di.notna().all()
    assert patho_di.min() > normal_di.max(), (
        f"expected pathological min ({patho_di.min()!r}) > normal max ({normal_di.max()!r})"
    )
