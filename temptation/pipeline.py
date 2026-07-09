"""The single entry point both CLI and GUI call. Neither front-end may
reimplement any of the orchestration below (see IMPLEMENTATION_BLUEPRINT.md
Sec 4.4 for the automated check that enforces this on the GUI)."""

import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from . import masks as masks_mod
from . import segmentation as seg_mod
from . import metrics_axon
from . import metrics_mito
from . import metrics_spatial
from .config import SegmentationConfig
from .schema import AXON_COLUMNS_V2, IMAGE_COLUMNS_V2, conform


def analyze_image_legacy(
    tem: np.ndarray,
    mask: np.ndarray,
    pixel_length_um: float,
    seg_cfg: SegmentationConfig,
):
    """Verbatim-equivalent port of the original measure_image() (Phase 1).

    Returns (df_axons, df_image, labels_ws, resolved_mode), exactly
    matching the legacy 4-tuple and legacy column set/order. `tem` is
    accepted but unused (matches the original -- see AUDIT.md / F8).
    """
    px2 = pixel_length_um ** 2

    # Stage 1: resolve mode (before compartments are built, matching the
    # original call order).
    resolved_mode = seg_mod.resolve_mode(mask, seg_cfg)

    # Stage 2: binary compartments & cleanup.
    comps = masks_mod.build_compartments(mask, seg_cfg)

    # Stage 3: axon seeds.
    axon_lab = seg_mod.label_axons(comps, seg_cfg)

    # Stage 4: watershed.
    labels_ws = seg_mod.run_watershed(axon_lab, comps, resolved_mode, seg_cfg)

    # Stage 4b: nearest-axon map for detached myelin, if requested.
    if seg_cfg.assign_detached_myelin == "nearest":
        nearest_axon_map_global = seg_mod.nearest_axon_map(axon_lab)
    else:
        nearest_axon_map_global = None

    # Stage 4c: global mito-to-axon assignment (Phase 3a, fixes F8 --
    # see segmentation.assign_mito_to_axons). "legacy" skips this and the
    # per-axon loop below falls back to the original per-fiber-crop
    # intersection instead.
    n_mito_assigned = np.nan
    n_mito_unassigned = np.nan
    fiber_to_mito_regions = {}
    if seg_cfg.mito_assignment != "legacy":
        mito_lab_global, mito_assignment_map = seg_mod.assign_mito_to_axons(
            comps.mito_in_axon, labels_ws, method=seg_cfg.mito_assignment,
        )
        mito_regions_by_label = {r.label: r for r in regionprops(mito_lab_global)}
        for mito_label, fiber_label in mito_assignment_map.items():
            fiber_to_mito_regions.setdefault(fiber_label, []).append(mito_regions_by_label[mito_label])
        n_mito_assigned = len(mito_assignment_map)
        n_mito_unassigned = int(mito_lab_global.max()) - n_mito_assigned

    # Stage 5: per-axon metric extraction.
    H, W = mask.shape[:2]
    fov_area_um2 = H * W * px2

    props = regionprops(labels_ws)
    rows = []

    for p in props:
        fiber_id = p.label
        min_row, min_col, max_row, max_col = p.bbox
        fiber_crop = p.image

        axon_crop = comps.axoplasm[min_row:max_row, min_col:max_col]
        axon_mask_local = fiber_crop & axon_crop

        axon_area_px = int(axon_mask_local.sum())
        if axon_area_px == 0:
            continue

        if seg_cfg.assign_detached_myelin == "nearest":
            myelin_area_px = int((comps.myelin & (nearest_axon_map_global == fiber_id)).sum())
            fiber_area_px = axon_area_px + myelin_area_px
        else:
            fiber_area_px = p.area
            myelin_area_px = fiber_area_px - axon_area_px

        geom = metrics_axon.axon_geometry(
            axon_area_px=axon_area_px,
            fiber_area_px=fiber_area_px,
            myelin_area_px=myelin_area_px,
            resolved_mode=resolved_mode,
            min_myelin_area_px=seg_cfg.min_myelin_area_px,
            pixel_length_um=pixel_length_um,
        )

        axon_labeled_local = label(axon_mask_local)
        axon_rp = regionprops(axon_labeled_local)
        if not axon_rp:
            continue
        axon_props_local = max(axon_rp, key=lambda r: r.area)

        shape = metrics_axon.shape_metrics(axon_props_local, axon_area_px, pixel_length_um)

        cy_local, cx_local = axon_props_local.centroid
        cy_global = min_row + cy_local
        cx_global = min_col + cx_local
        cx_um = cx_global * pixel_length_um
        cy_um = cy_global * pixel_length_um

        if seg_cfg.mito_assignment == "legacy":
            mito_crop = comps.mito[min_row:max_row, min_col:max_col]
            mito_in_axon_crop = mito_crop & axon_mask_local
            mito = metrics_mito.mito_metrics_for_axon(
                mito_in_axon_crop=mito_in_axon_crop,
                axon_area_um2=geom["axon_area_um2"],
                cx_local=cx_local,
                cy_local=cy_local,
                pixel_length_um=pixel_length_um,
            )
        else:
            mito = metrics_mito.mito_metrics_from_regions(
                mito_regions=fiber_to_mito_regions.get(fiber_id, []),
                axon_area_um2=geom["axon_area_um2"],
                axon_cx_px=cx_global,
                axon_cy_px=cy_global,
                pixel_length_um=pixel_length_um,
            )

        rows.append({
            "axon_id": fiber_id,
            "axon_area_um2": geom["axon_area_um2"],
            "myelin_area_um2": geom["myelin_area_um2"],
            "fiber_area_um2": geom["fiber_area_um2"],
            "d_inner_um": geom["d_inner_um"],
            "d_outer_um": geom["d_outer_um"],
            "g_ratio": geom["g_ratio"],
            "myelin_thickness_um": geom["myelin_thickness_um"],
            "perimeter_um": shape["perimeter_um"],
            "eccentricity": shape["eccentricity"],
            "solidity": shape["solidity"],
            "circularity": shape["circularity"],
            "convexity": shape["convexity"],
            "axon_vol_fraction": geom["axon_vol_fraction"],
            "myelin_vol_fraction": geom["myelin_vol_fraction"],
            "centroid_x_um": cx_um,
            "centroid_y_um": cy_um,
            "centroid_x_px": cx_global,
            "centroid_y_px": cy_global,
            "mito_count": mito["mito_count"],
            "mito_area_um2": mito["mito_area_um2"],
            "mito_density_per_um2": mito["mito_density_per_um2"],
            "mito_mean_circularity": mito["mito_mean_circularity"],
            "mito_mean_form_factor": mito["mito_mean_form_factor"],
            "mito_mean_feret_um": mito["mito_mean_feret_um"],
            "mito_std_area_um2": mito["mito_std_area_um2"],
            "mito_max_area_um2": mito["mito_max_area_um2"],
            "mito_cv_area": mito["mito_cv_area"],
            "mito_area_skewness": mito["mito_area_skewness"],
            "mito_mean_dist_centroid_um": mito["mito_mean_dist_centroid_um"],
            "mito_std_dist_centroid_um": mito["mito_std_dist_centroid_um"],
            "axon_circularity_crofton": shape["axon_circularity_crofton"],
            "axon_shape_irregularity": shape["axon_shape_irregularity"],
            "axon_shape_irregularity_crofton": shape["axon_shape_irregularity_crofton"],
        })

    df_axons = pd.DataFrame(rows)

    # Stage 6: nearest-neighbour distances.
    if not df_axons.empty:
        if len(df_axons) >= 2:
            centroids = df_axons[["centroid_x_um", "centroid_y_um"]].to_numpy()
            df_axons["nearest_neighbor_um"] = metrics_spatial.axon_nearest_neighbors(centroids)
        else:
            df_axons["nearest_neighbor_um"] = np.nan

    # Stage 7: image-level summary.
    if not df_axons.empty:
        n_fibers = len(df_axons)
        fiber_density_per_mm2 = (n_fibers / fov_area_um2) * 1e6

        summary = {
            "n_fibers": n_fibers,
            "fov_area_um2": fov_area_um2,
            "fiber_density_per_mm2": fiber_density_per_mm2,
            "mean_axon_area_um2": df_axons["axon_area_um2"].mean(),
            "mean_circularity": df_axons["circularity"].mean(),
            "mean_convexity": df_axons["convexity"].mean(),
            "mean_avf": df_axons["axon_vol_fraction"].mean(),
            "mode": resolved_mode,
        }
        if resolved_mode == "normal":
            summary["mean_g_ratio"] = df_axons["g_ratio"].mean()
            summary["mean_myelin_thickness_um"] = df_axons["myelin_thickness_um"].mean()
            summary["mean_mvf"] = df_axons["myelin_vol_fraction"].mean()
        else:
            summary["mean_g_ratio"] = np.nan
            summary["mean_myelin_thickness_um"] = np.nan
            summary["mean_mvf"] = np.nan

        # F2 (IMPLEMENTATION_BLUEPRINT.md Sec 0): the legacy computation
        # compares mito against axon_only (raw, pre-hole-handling), which
        # never overlaps mito by construction -- mito_outside_frac was
        # therefore *always* exactly 1.0, regardless of mode. Preserved
        # exactly under "legacy" for regression; under "fill" it correctly
        # compares against the already-fixed axoplasm (which now includes
        # every mitochondrion actually attributed to an axon), so the
        # metric becomes meaningful again.
        if seg_cfg.mito_hole_handling == "legacy":
            mito_outside_reference = comps.axon_only
        else:
            mito_outside_reference = comps.axoplasm
        mito_outside = comps.mito & ~(mito_outside_reference | comps.myelin)
        mito_outside_area_um2 = float(mito_outside.sum()) * px2
        total_mito_px = int(comps.mito.sum())
        mito_outside_frac = float(mito_outside.sum()) / max(total_mito_px, 1)
        summary["mito_outside_area_um2"] = mito_outside_area_um2
        summary["mito_outside_frac"] = mito_outside_frac

        total_myelin_px = int((mask == seg_cfg.myelin_val).sum())
        total_myelin_um2 = total_myelin_px * px2
        myelin_area_fraction_of_fov = total_myelin_um2 / fov_area_um2 if fov_area_um2 > 0 else 0.0
        myelin_lab = label(mask == seg_cfg.myelin_val)
        myelin_component_count = int(myelin_lab.max())
        summary["total_myelin_area_um2"] = total_myelin_um2
        summary["myelin_area_fraction_of_fov"] = myelin_area_fraction_of_fov
        summary["myelin_component_count"] = myelin_component_count

        # Phase 3a diagnostic: NaN under mito_assignment="legacy" (no
        # global assignment is computed in that mode -- see Stage 4c).
        summary["n_mito_assigned"] = n_mito_assigned
        summary["n_mito_unassigned"] = n_mito_unassigned

        df_image = pd.DataFrame([summary])
    else:
        df_image = pd.DataFrame()

    if not df_axons.empty:
        df_axons = conform(df_axons, [c for c in AXON_COLUMNS_V2 if c not in ("image_id", "mode")])
    if not df_image.empty:
        df_image = conform(df_image, [c for c in IMAGE_COLUMNS_V2 if c != "image_id"])

    return df_axons, df_image, labels_ws, resolved_mode
