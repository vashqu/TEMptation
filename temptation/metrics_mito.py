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
