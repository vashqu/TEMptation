"""Inter-axon spatial arrangement (image-scoped, not per-axon-region)."""

import numpy as np
from scipy.spatial import cKDTree


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
