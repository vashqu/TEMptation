"""Spatial arrangement metrics: inter-axon (whole image) and
intra-axon/mitochondrial (Phase 3b, per axon)."""

import numpy as np
from scipy.spatial import cKDTree

from .mathutils import equivalent_radius, safe_divide


def axon_nearest_neighbors(centroids_um: np.ndarray) -> np.ndarray:
    """Distance from each axon centroid to its nearest other axon centroid.

    Verbatim port of measure_nerve.py:435-442. Caller is responsible for
    only invoking this when len(centroids_um) >= 2 (matching the original
    guard); this function always uses the k=2 cKDTree query.
    """
    centroids_um = np.asarray(centroids_um)
    tree = cKDTree(centroids_um)
    dists, _ = tree.query(centroids_um, k=2)
    return dists[:, 1]


def mito_peripheralization_index(mean_dist_centroid_um: float, axon_area_um2: float) -> float:
    """CLAUDE.md Sec 5.4.1: mean distance of mitochondria from the axon
    center, normalized by an equivalent axon radius. ~0 = central,
    ~1 = at the boundary; values above 1 can occur for irregular axons and
    should be flagged (Phase 5 QC), not discarded.

    A pure derived ratio of two already-computed scalars -- does not need
    raw mitochondrion data, so it works identically regardless of which
    mito_assignment mode produced mean_dist_centroid_um.
    """
    return safe_divide(mean_dist_centroid_um, equivalent_radius(axon_area_um2))


def mito_spatial_clustering(
    mito_regions: list,
    axon_area_um2: float,
    pixel_length_um: float,
    min_count_for_clustering: int = 3,
) -> dict:
    """CLAUDE.md Sec 5.4.2: nearest-neighbor clustering among mitochondria
    within one axon.

    mito_mean_nn_distance_um needs >=2 mitochondria (a NN distance is
    undefined for 0 or 1 points) and is reported whenever available.
    mito_clustering_index additionally needs >=min_count_for_clustering
    (default 3, per CLAUDE.md's explicit caution that the Poisson-process
    approximation is unstable below that) and is NaN otherwise even when
    the NN distance itself is computable.

    Index ~1: approximately random (CSR). >1: observed NN distance is
    smaller than expected under CSR -- clustering. <1: dispersion.
    """
    mito_count = len(mito_regions)

    if mito_count < 2:
        return {"mito_mean_nn_distance_um": np.nan, "mito_clustering_index": np.nan}

    centroids_px = np.array([mr.centroid for mr in mito_regions])
    centroids_um = centroids_px * pixel_length_um
    tree = cKDTree(centroids_um)
    dists, _ = tree.query(centroids_um, k=2)
    observed_mean_nn_um = float(np.mean(dists[:, 1]))

    if mito_count < min_count_for_clustering:
        return {"mito_mean_nn_distance_um": observed_mean_nn_um, "mito_clustering_index": np.nan}

    lam = safe_divide(mito_count, axon_area_um2)
    if np.isnan(lam) or lam <= 0:
        clustering_index = np.nan
    else:
        expected_nn_um = 1.0 / (2.0 * np.sqrt(lam))
        clustering_index = safe_divide(expected_nn_um, observed_mean_nn_um)

    return {
        "mito_mean_nn_distance_um": observed_mean_nn_um,
        "mito_clustering_index": clustering_index,
    }
