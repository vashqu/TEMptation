"""Per-axon mitochondrial aggregate metrics.

Two parallel functions, deliberately not unified:

- `mito_metrics_for_axon` -- the Phase 1 legacy path. Re-labels a
  per-fiber crop each call (matching the original measure_image()
  exactly) and operates in crop-local pixel coordinates. Regression-
  tested against the Phase 0/2b golden baselines; never modify it.
- `mito_metrics_from_regions` -- Phase 3. Takes pre-labeled GLOBAL
  regionprops (from segmentation.assign_mito_to_axons, called once per
  image rather than once per axon) and operates in global pixel
  coordinates. This is what fixes F8 (a mitochondrion straddling a
  watershed boundary was being fragmented and double-counted by the
  crop-based approach).
"""

import numpy as np
from scipy.stats import skew as scipy_skew
from skimage.measure import label, regionprops

from .mathutils import safe_aspect_ratio, safe_divide


def mito_metrics_for_axon(
    mito_in_axon_crop: np.ndarray,
    axon_area_um2: float,
    cx_local: float,
    cy_local: float,
    pixel_length_um: float,
) -> dict:
    """Verbatim port of measure_nerve.py:302-385."""
    px2 = pixel_length_um ** 2

    mito_area_px_loc = int(mito_in_axon_crop.sum())
    mito_area_um2 = mito_area_px_loc * px2
    mito_lab = label(mito_in_axon_crop, connectivity=2)
    mito_regions = regionprops(mito_lab)
    mito_count = len(mito_regions)
    mito_density = mito_count / axon_area_um2 if axon_area_um2 > 0 else np.nan

    if mito_count > 0:
        mito_circularities = []
        mito_form_factors = []
        mito_ferets_um = []
        mito_areas_um2 = []
        mito_cx_list = []
        mito_cy_list = []

        for mr in mito_regions:
            m_area = mr.area
            m_perim = mr.perimeter
            m_areas_um2 = m_area * px2
            mito_areas_um2.append(m_areas_um2)

            if m_perim > 0:
                m_circ = (4.0 * np.pi * m_area) / (m_perim ** 2)
            else:
                m_circ = np.nan
            mito_circularities.append(m_circ)

            if m_perim > 0 and m_area > 0:
                m_ff = (m_perim ** 2) / (4.0 * np.pi * m_area)
            else:
                m_ff = np.nan
            mito_form_factors.append(m_ff)

            m_feret_um = mr.feret_diameter_max * pixel_length_um
            mito_ferets_um.append(m_feret_um)

            mito_cy_list.append(mr.centroid[0])
            mito_cx_list.append(mr.centroid[1])

        mito_areas_arr = np.array(mito_areas_um2)

        mito_mean_circularity = float(np.nanmean(mito_circularities))
        mito_mean_form_factor = float(np.nanmean(mito_form_factors))
        mito_mean_feret_um = float(np.nanmean(mito_ferets_um))
        mito_std_area_um2 = float(np.std(mito_areas_arr)) if mito_count > 1 else 0.0
        mito_max_area_um2 = float(np.max(mito_areas_arr))
        mito_cv_area = (
            mito_std_area_um2 / float(np.mean(mito_areas_arr))
            if mito_count > 1 and np.mean(mito_areas_arr) > 0
            else np.nan
        )
        mito_area_skewness = (
            float(scipy_skew(mito_areas_arr))
            if mito_count > 2
            else np.nan
        )

        mito_cx_arr = np.array(mito_cx_list)
        mito_cy_arr = np.array(mito_cy_list)
        dists_px = np.sqrt((mito_cx_arr - cx_local) ** 2 +
                            (mito_cy_arr - cy_local) ** 2)
        dists_um = dists_px * pixel_length_um
        mito_mean_dist_centroid_um = float(np.mean(dists_um))
        mito_std_dist_centroid_um = float(np.std(dists_um)) if mito_count > 1 else 0.0

    else:
        mito_mean_circularity = np.nan
        mito_mean_form_factor = np.nan
        mito_mean_feret_um = np.nan
        mito_std_area_um2 = np.nan
        mito_max_area_um2 = np.nan
        mito_cv_area = np.nan
        mito_area_skewness = np.nan
        mito_mean_dist_centroid_um = np.nan
        mito_std_dist_centroid_um = np.nan

    return {
        "mito_count": mito_count,
        "mito_area_um2": mito_area_um2,
        "mito_density_per_um2": mito_density,
        "mito_mean_circularity": mito_mean_circularity,
        "mito_mean_form_factor": mito_mean_form_factor,
        "mito_mean_feret_um": mito_mean_feret_um,
        "mito_std_area_um2": mito_std_area_um2,
        "mito_max_area_um2": mito_max_area_um2,
        "mito_cv_area": mito_cv_area,
        "mito_area_skewness": mito_area_skewness,
        "mito_mean_dist_centroid_um": mito_mean_dist_centroid_um,
        "mito_std_dist_centroid_um": mito_std_dist_centroid_um,
    }


def mito_metrics_from_regions(
    mito_regions: list,
    axon_area_um2: float,
    axon_cx_px: float,
    axon_cy_px: float,
    pixel_length_um: float,
) -> dict:
    """Phase 3: same aggregate output as mito_metrics_for_axon, but from a
    pre-selected list of GLOBAL regionprops (produced by
    segmentation.assign_mito_to_axons, one full/unfragmented mitochondrion
    per entry) instead of re-labeling a per-fiber crop. Kept as a fully
    separate function -- not a refactor of mito_metrics_for_axon -- so the
    regression-tested legacy path (IMPLEMENTATION_BLUEPRINT.md Sec 9
    Phase 1/2b) is never touched by this addition.

    `mito_regions[i].centroid` and `axon_cx_px`/`axon_cy_px` must all be in
    the SAME (global image) pixel coordinate frame -- unlike
    mito_metrics_for_axon, which operates in per-fiber-crop-local
    coordinates throughout.
    """
    px2 = pixel_length_um ** 2
    mito_count = len(mito_regions)
    mito_area_um2 = sum(mr.area for mr in mito_regions) * px2
    mito_density = mito_count / axon_area_um2 if axon_area_um2 > 0 else np.nan

    if mito_count > 0:
        mito_circularities = []
        mito_form_factors = []
        mito_ferets_um = []
        mito_areas_um2 = []
        mito_cx_list = []
        mito_cy_list = []

        for mr in mito_regions:
            m_area = mr.area
            m_perim = mr.perimeter
            m_areas_um2 = m_area * px2
            mito_areas_um2.append(m_areas_um2)

            if m_perim > 0:
                m_circ = (4.0 * np.pi * m_area) / (m_perim ** 2)
            else:
                m_circ = np.nan
            mito_circularities.append(m_circ)

            if m_perim > 0 and m_area > 0:
                m_ff = (m_perim ** 2) / (4.0 * np.pi * m_area)
            else:
                m_ff = np.nan
            mito_form_factors.append(m_ff)

            m_feret_um = mr.feret_diameter_max * pixel_length_um
            mito_ferets_um.append(m_feret_um)

            mito_cy_list.append(mr.centroid[0])
            mito_cx_list.append(mr.centroid[1])

        mito_areas_arr = np.array(mito_areas_um2)

        mito_mean_circularity = float(np.nanmean(mito_circularities))
        mito_mean_form_factor = float(np.nanmean(mito_form_factors))
        mito_mean_feret_um = float(np.nanmean(mito_ferets_um))
        mito_std_area_um2 = float(np.std(mito_areas_arr)) if mito_count > 1 else 0.0
        mito_max_area_um2 = float(np.max(mito_areas_arr))
        mito_cv_area = (
            mito_std_area_um2 / float(np.mean(mito_areas_arr))
            if mito_count > 1 and np.mean(mito_areas_arr) > 0
            else np.nan
        )
        mito_area_skewness = (
            float(scipy_skew(mito_areas_arr))
            if mito_count > 2
            else np.nan
        )

        mito_cx_arr = np.array(mito_cx_list)
        mito_cy_arr = np.array(mito_cy_list)
        dists_px = np.sqrt((mito_cx_arr - axon_cx_px) ** 2 +
                            (mito_cy_arr - axon_cy_px) ** 2)
        dists_um = dists_px * pixel_length_um
        mito_mean_dist_centroid_um = float(np.mean(dists_um))
        mito_std_dist_centroid_um = float(np.std(dists_um)) if mito_count > 1 else 0.0

    else:
        mito_mean_circularity = np.nan
        mito_mean_form_factor = np.nan
        mito_mean_feret_um = np.nan
        mito_std_area_um2 = np.nan
        mito_max_area_um2 = np.nan
        mito_cv_area = np.nan
        mito_area_skewness = np.nan
        mito_mean_dist_centroid_um = np.nan
        mito_std_dist_centroid_um = np.nan

    return {
        "mito_count": mito_count,
        "mito_area_um2": mito_area_um2,
        "mito_density_per_um2": mito_density,
        "mito_mean_circularity": mito_mean_circularity,
        "mito_mean_form_factor": mito_mean_form_factor,
        "mito_mean_feret_um": mito_mean_feret_um,
        "mito_std_area_um2": mito_std_area_um2,
        "mito_max_area_um2": mito_max_area_um2,
        "mito_cv_area": mito_cv_area,
        "mito_area_skewness": mito_area_skewness,
        "mito_mean_dist_centroid_um": mito_mean_dist_centroid_um,
        "mito_std_dist_centroid_um": mito_std_dist_centroid_um,
    }


def mito_burden_and_shape_metrics(
    mito_regions: list,
    axon_area_um2: float,
    fiber_area_um2: float,
    myelin_area_um2: float,
    pixel_length_um: float,
) -> dict:
    """Phase 3b per-axon burden and shape metrics (occupancy, size
    distribution, fragmentation, load relative to fiber/myelin, and
    aspect-ratio/solidity/eccentricity aggregates).

    Deliberately takes a raw `mito_regions` list rather than a mask or a
    crop -- it doesn't care whether the regions came from the legacy
    per-fiber-crop labeling or the Phase 3a global assignment, so the same
    function serves both `mito_assignment` modes. New metrics use
    ddof=1 (see mathutils.py), unlike the legacy std/cv columns above,
    which use ddof=0.

    Empty-mito contract (CLAUDE.md Sec 12.3): count=0 -> total area 0,
    occupancy/load ratios 0 (well-defined: 0/positive-denominator),
    fragmentation and all shape aggregates NaN (undefined, not 0/0).
    """
    px2 = pixel_length_um ** 2
    mito_count = len(mito_regions)
    areas_um2 = np.array([mr.area * px2 for mr in mito_regions], dtype=float)
    mito_total_area_um2 = float(areas_um2.sum()) if mito_count > 0 else 0.0

    mito_occupancy_ratio = safe_divide(mito_total_area_um2, axon_area_um2)
    normalized_mito_load = safe_divide(mito_total_area_um2, fiber_area_um2)
    mito_per_myelin = safe_divide(mito_total_area_um2, myelin_area_um2)
    mito_fragmentation_index = (
        safe_divide(mito_count, mito_total_area_um2) if mito_count > 0 else np.nan
    )

    if mito_count > 0:
        mito_mean_area_um2 = float(np.mean(areas_um2))
        mito_median_area_um2 = float(np.median(areas_um2))
        q75, q25 = np.percentile(areas_um2, [75, 25])
        mito_area_iqr = float(q75 - q25)

        aspect_ratios = np.array([safe_aspect_ratio(mr) for mr in mito_regions], dtype=float)
        solidities = np.array([mr.solidity for mr in mito_regions], dtype=float)
        eccentricities = np.array([mr.eccentricity for mr in mito_regions], dtype=float)

        with np.errstate(invalid="ignore"):
            mito_mean_aspect_ratio = float(np.nanmean(aspect_ratios))
            mito_std_aspect_ratio = (
                float(np.nanstd(aspect_ratios, ddof=1)) if mito_count > 1 else np.nan
            )
            mito_mean_solidity = float(np.nanmean(solidities))
            mito_std_solidity = (
                float(np.nanstd(solidities, ddof=1)) if mito_count > 1 else np.nan
            )
            mito_mean_eccentricity = float(np.nanmean(eccentricities))
    else:
        mito_mean_area_um2 = np.nan
        mito_median_area_um2 = np.nan
        mito_area_iqr = np.nan
        mito_mean_aspect_ratio = np.nan
        mito_std_aspect_ratio = np.nan
        mito_mean_solidity = np.nan
        mito_std_solidity = np.nan
        mito_mean_eccentricity = np.nan

    return {
        "mito_total_area_um2": mito_total_area_um2,
        "mito_occupancy_ratio": mito_occupancy_ratio,
        "mito_mean_area_um2": mito_mean_area_um2,
        "mito_median_area_um2": mito_median_area_um2,
        "mito_area_iqr": mito_area_iqr,
        "mito_fragmentation_index": mito_fragmentation_index,
        "normalized_mito_load": normalized_mito_load,
        "mito_per_myelin": mito_per_myelin,
        "mito_mean_aspect_ratio": mito_mean_aspect_ratio,
        "mito_std_aspect_ratio": mito_std_aspect_ratio,
        "mito_mean_solidity": mito_mean_solidity,
        "mito_std_solidity": mito_std_solidity,
        "mito_mean_eccentricity": mito_mean_eccentricity,
    }
