"""The single entry point both CLI and GUI call. Neither front-end may
reimplement any of the orchestration below (see IMPLEMENTATION_BLUEPRINT.md
Sec 4.4 for the automated check that enforces this on the GUI)."""

import traceback
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from . import dataio
from . import masks as masks_mod
from . import qc as qc_mod
from . import segmentation as seg_mod
from . import metrics_axon
from . import metrics_mito
from . import metrics_spatial
from . import summaries
from .config import QCThresholds, SegmentationConfig
from .mathutils import safe_aspect_ratio
from .schema import AXON_COLUMNS_V2, IMAGE_COLUMNS_V2, MITO_COLUMNS, conform, resolve_schema_version


def analyze_image_legacy(
    tem: np.ndarray,
    mask: np.ndarray,
    pixel_length_um: float,
    seg_cfg: SegmentationConfig,
    reference_myelin_fraction: float = None,
    qc_thresholds: QCThresholds = None,
    exclude_qc_failed: bool = False,
):
    """Started as a verbatim-equivalent port of the original measure_image()
    (Phase 1); the axon/image column *content* for the legacy fields is
    still byte-exact (see tests/golden/), but the return signature grew a
    5th element in Phase 3b.

    Returns (df_axons, df_image, labels_ws, resolved_mode, df_mito).
    df_mito is empty when mito_assignment="legacy" (no global per-
    mitochondrion table is meaningful in that mode -- see
    docs/phase3a_impact.md). compat.measure_image() unpacks this 5-tuple
    and returns only the first 4, preserving the legacy signature for the
    GUI. `tem` is accepted but unused (matches the original -- see
    AUDIT.md / F8).

    `reference_myelin_fraction` (Phase 4) feeds image_demyelination_index
    (see summaries.demyelination_index). Deliberately no default: this is
    a biological calibration choice IMPLEMENTATION_BLUEPRINT.md Sec 11 A4
    flagged as not derivable from the code or data, and the user was
    asked explicitly and chose no default over a self-calibrating one --
    image_demyelination_index is NaN unless this is supplied.

    `qc_thresholds`/`exclude_qc_failed` (Phase 5): qc_* flags are always
    computed (defaults to QCThresholds() if not given -- flags are cheap
    and CLAUDE.md Sec 7 wants them on by default). exclude_qc_failed
    defaults to False, matching CLAUDE.md Sec 8's default behavior:
    compute everything, flag everything, exclude nothing unless the
    caller opts in. Excluded rows are marked, never dropped -- axons.csv
    row count is invariant to this flag (see tests/test_qc.py).
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
    mito_rows = []  # Phase 3b: one row per real mitochondrion, only
    # populated when mito_assignment != "legacy" (a "one row per
    # mitochondrion" table is only meaningful once nothing can fragment
    # into two rows -- see docs/phase3a_impact.md).

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
            # Phase 3b's burden/shape/spatial functions take a raw region
            # list regardless of assignment mode (see metrics_mito.py's
            # module docstring), so re-derive it here for the legacy path
            # without touching the regression-tested mito_metrics_for_axon.
            mito_regions_this_axon = regionprops(label(mito_in_axon_crop, connectivity=2))
        else:
            mito_regions_this_axon = fiber_to_mito_regions.get(fiber_id, [])
            mito = metrics_mito.mito_metrics_from_regions(
                mito_regions=mito_regions_this_axon,
                axon_area_um2=geom["axon_area_um2"],
                axon_cx_px=cx_global,
                axon_cy_px=cy_global,
                pixel_length_um=pixel_length_um,
            )

        burden = metrics_mito.mito_burden_and_shape_metrics(
            mito_regions=mito_regions_this_axon,
            axon_area_um2=geom["axon_area_um2"],
            fiber_area_um2=geom["fiber_area_um2"],
            myelin_area_um2=geom["myelin_area_um2"],
            pixel_length_um=pixel_length_um,
        )
        peripheralization = metrics_spatial.mito_peripheralization_index(
            mito["mito_mean_dist_centroid_um"], geom["axon_area_um2"],
        )
        clustering = metrics_spatial.mito_spatial_clustering(
            mito_regions_this_axon, geom["axon_area_um2"], pixel_length_um,
            min_count_for_clustering=3,
        )

        if seg_cfg.mito_assignment != "legacy" and mito_regions_this_axon:
            for mr in mito_regions_this_axon:
                mr_cy, mr_cx = mr.centroid
                dist_um = float(np.hypot(mr_cx - cx_global, mr_cy - cy_global)) * pixel_length_um
                mito_rows.append({
                    "mito_id": mr.label,
                    "parent_axon_id": fiber_id,
                    "mito_area_um2": mr.area * pixel_length_um ** 2,
                    "mito_perimeter_um": mr.perimeter * pixel_length_um,
                    "mito_aspect_ratio": safe_aspect_ratio(mr),
                    "mito_solidity": mr.solidity,
                    "mito_eccentricity": mr.eccentricity,
                    "mito_circularity": (
                        (4.0 * np.pi * mr.area) / (mr.perimeter ** 2) if mr.perimeter > 0 else np.nan
                    ),
                    "mito_feret_um": mr.feret_diameter_max * pixel_length_um,
                    "mito_centroid_x_px": mr_cx,
                    "mito_centroid_y_px": mr_cy,
                    "mito_dist_to_axon_center_um": dist_um,
                })

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
            "mito_total_area_um2": burden["mito_total_area_um2"],
            "mito_occupancy_ratio": burden["mito_occupancy_ratio"],
            "mito_mean_area_um2": burden["mito_mean_area_um2"],
            "mito_median_area_um2": burden["mito_median_area_um2"],
            "mito_area_iqr": burden["mito_area_iqr"],
            "mito_fragmentation_index": burden["mito_fragmentation_index"],
            "normalized_mito_load": burden["normalized_mito_load"],
            "mito_per_myelin": burden["mito_per_myelin"],
            "mito_mean_aspect_ratio": burden["mito_mean_aspect_ratio"],
            "mito_std_aspect_ratio": burden["mito_std_aspect_ratio"],
            "mito_mean_solidity": burden["mito_mean_solidity"],
            "mito_std_solidity": burden["mito_std_solidity"],
            "mito_mean_eccentricity": burden["mito_mean_eccentricity"],
            "mito_peripheralization_index": peripheralization,
            "mito_mean_nn_distance_um": clustering["mito_mean_nn_distance_um"],
            "mito_clustering_index": clustering["mito_clustering_index"],
        })

    df_axons = pd.DataFrame(rows)
    df_mito = pd.DataFrame(mito_rows)
    # Phase 6 (fixes the unstable-schema half of F8): conform
    # unconditionally, even when 0 rows, so a 0-mitochondrion image still
    # yields the full MITO_COLUMNS header rather than an empty (0-column)
    # frame -- verified this matters: pd.DataFrame().reindex(columns=...)
    # correctly produces a 0-row, fully-columned frame, and concatenating
    # multiple such frames across a batch preserves that column set even
    # if every single image in the batch has zero mitochondria.
    df_mito = conform(df_mito, MITO_COLUMNS)

    # Stage 6: nearest-neighbour distances.
    if not df_axons.empty:
        if len(df_axons) >= 2:
            centroids = df_axons[["centroid_x_um", "centroid_y_um"]].to_numpy()
            df_axons["nearest_neighbor_um"] = metrics_spatial.axon_nearest_neighbors(centroids)
        else:
            df_axons["nearest_neighbor_um"] = np.nan

    # Stage 6b (Phase 5): QC flags (always computed) and exclusion marking
    # (opt-in, never drops rows -- CLAUDE.md Sec 6/8).
    effective_qc_thresholds = qc_thresholds if qc_thresholds is not None else QCThresholds()
    df_axons = qc_mod.compute_qc_flags(df_axons, effective_qc_thresholds)
    df_axons = qc_mod.apply_exclusions(df_axons, enabled=exclude_qc_failed)

    # Stage 7: image-level summary.
    if not df_axons.empty:
        n_fibers = len(df_axons)
        fiber_density_per_mm2 = (n_fibers / fov_area_um2) * 1e6

        valid_mask = ~df_axons["excluded_from_analysis"]
        df_axons_valid = df_axons[valid_mask]
        n_valid = int(valid_mask.sum())
        n_excluded = n_fibers - n_valid

        summary = {
            "n_fibers": n_fibers,
            "fov_area_um2": fov_area_um2,
            "fiber_density_per_mm2": fiber_density_per_mm2,
            "mode": resolved_mode,
        }
        # Axon-row-dependent aggregates: main columns reflect the QC-valid
        # subset (identical to "_all" when nothing is excluded, which is
        # why this reproduces the legacy mean_* values exactly by default
        # -- see tests/test_regression_golden.py). "_all" always present
        # for schema stability and to preserve the unfiltered view
        # (CLAUDE.md Sec 6.4).
        summary.update(summaries.summarize_axon_aggregates(df_axons_valid, resolved_mode))
        summary_all = summaries.summarize_axon_aggregates(df_axons, resolved_mode)
        summary.update({f"{k}_all": v for k, v in summary_all.items()})

        summary["image_n_axons_total"] = n_fibers
        summary["image_n_axons_valid"] = n_valid
        summary["image_n_axons_excluded"] = n_excluded
        flag_cols = [c for c in qc_mod.QC_FLAG_COLUMNS if c in df_axons.columns]
        any_flag = df_axons[flag_cols].any(axis=1) if flag_cols else pd.Series(False, index=df_axons.index)
        summary["image_percent_flagged_axons"] = (
            100.0 * any_flag.mean() if n_fibers > 0 else np.nan
        )
        summary["qc_high_exclusion_rate"] = (
            bool(n_excluded / n_fibers > 0.3) if n_fibers > 0 else False
        )

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

        # Phase 4: demyelination index (mask-derived, not axon-row-
        # dependent, so computed once regardless of exclusion).
        summary["image_demyelination_index"] = summaries.demyelination_index(
            myelin_area_fraction_of_fov, reference_myelin_fraction,
        )

        # Phase 5: mask-quality QC diagnostic (F6) -- wires up
        # masks.mask_sanity, unused since Phase 2.
        mask_qc = masks_mod.mask_sanity(mask, seg_cfg)
        summary["image_noncanonical_mask_frac"] = mask_qc["noncanonical_frac"]
        summary["qc_mask_noncanonical"] = bool(mask_qc["noncanonical_frac"] > 0.001)

        df_image = pd.DataFrame([summary])
    else:
        # Legacy behavior preserved: an image with zero detected axons
        # gets zero image_summary.csv rows (not a placeholder row of
        # NaNs) -- matching the original measure_image()'s "if not
        # df_axons.empty" gate. Phase 6 only fixes the *columns* of that
        # empty result, not whether a row is emitted at all.
        df_image = pd.DataFrame()

    # Phase 6: conform unconditionally (see the df_mito comment above for
    # why this is safe and necessary even at 0 rows).
    df_axons = conform(df_axons, [c for c in AXON_COLUMNS_V2 if c not in ("image_id", "mode")])
    df_image = conform(df_image, [c for c in IMAGE_COLUMNS_V2 if c != "image_id"])

    return df_axons, df_image, labels_ws, resolved_mode, df_mito


@dataclass
class DatasetResult:
    df_axons: pd.DataFrame = field(default_factory=pd.DataFrame)
    df_image: pd.DataFrame = field(default_factory=pd.DataFrame)
    df_mito: pd.DataFrame = field(default_factory=pd.DataFrame)
    df_qc_report: pd.DataFrame = field(default_factory=pd.DataFrame)
    errors: list = field(default_factory=list)  # [(pair_id, message, traceback_str), ...]


def analyze_dataset(
    pairs,
    pixel_length_um: float,
    seg_cfg: SegmentationConfig,
    qc_thresholds: QCThresholds = None,
    exclude_qc_failed: bool = False,
    reference_myelin_fraction: float = None,
    collect_mito: bool = False,
    collect_qc_report: bool = False,
    on_image_done=None,
    on_image_error=None,
    cancel_event=None,
) -> DatasetResult:
    """Batch orchestration shared by CLI and GUI (Phase 7 prep --
    IMPLEMENTATION_BLUEPRINT.md Sec 4.3's original "analyze_dataset"
    contract). For each pair: reads tem/mask, calls analyze_image_legacy,
    stamps identity/metadata columns (image_id, mode, group, image_path,
    mask_path, pixel_size_um, schema_version -- moved here from cli.py's
    former process_pair, since this is shared orchestration, not a
    CLI-specific concern), and accumulates.

    `on_image_done(index, total, pair, df_axons, df_image, labels_ws,
    resolved_mode, df_mito, tem, mask)` fires after each successfully
    processed image. `on_image_error(index, total, pair, exc, traceback_str)`
    fires on a per-image failure (caught, not fatal to the batch -- matches
    the original per-pair try/except in cli.py's main()). Neither
    printing nor plotting happens in this function itself -- cli.py's
    callback handles --plot and console output; the GUI's callback
    updates its own progress UI. This keeps analyze_dataset decoupled
    from what any particular caller wants to do with a per-image result.

    `cancel_event`, if given and set mid-run (anything with an
    `is_set()` method, typically threading.Event), stops processing
    further pairs; results already accumulated are kept, not discarded.

    `collect_mito`/`collect_qc_report` gate whether the (otherwise
    per-image, always computed) mito table / QC report get accumulated
    into the batch-level DatasetResult -- set False when the caller has
    no intention of writing/displaying them, to skip the accumulation
    work on large batches.
    """
    all_axons = []
    all_images = []
    all_mito = []
    all_qc_reports = []
    errors = []

    total = len(pairs)
    for i, pair in enumerate(pairs):
        if cancel_event is not None and cancel_event.is_set():
            break
        try:
            tem = dataio.read_image(pair["tem_path"])
            mask = dataio.read_mask(pair["mask_path"])

            df_axons, df_image, labels_ws, resolved_mode, df_mito = analyze_image_legacy(
                tem, mask, pixel_length_um, seg_cfg,
                reference_myelin_fraction=reference_myelin_fraction,
                qc_thresholds=qc_thresholds,
                exclude_qc_failed=exclude_qc_failed,
            )

            image_id = pair["id"]
            group = pair.get("group")
            for df in (df_axons, df_image):
                if "image_id" not in df.columns:
                    df.insert(0, "image_id", image_id)
                else:
                    df["image_id"] = image_id
                if "mode" not in df.columns:
                    df.insert(1, "mode", resolved_mode)
                else:
                    df["mode"] = resolved_mode
                df["group"] = group
                df["image_path"] = str(pair["tem_path"])
                df["mask_path"] = str(pair["mask_path"])
                df["pixel_size_um"] = pixel_length_um
                df["schema_version"] = resolve_schema_version(seg_cfg.mito_hole_handling)

            if collect_mito and not df_mito.empty:
                df_mito.insert(0, "image_id", image_id)
                df_mito["group"] = group

            df_qc_report = qc_mod.qc_report(df_axons) if collect_qc_report else pd.DataFrame()

            all_axons.append(df_axons)
            all_images.append(df_image)
            if collect_mito and not df_mito.empty:
                all_mito.append(df_mito)
            if not df_qc_report.empty:
                all_qc_reports.append(df_qc_report)

            if on_image_done is not None:
                on_image_done(i, total, pair, df_axons, df_image, labels_ws, resolved_mode, df_mito, tem, mask)

        except Exception as exc:
            tb_str = traceback.format_exc()
            errors.append((pair.get("id", "?"), str(exc), tb_str))
            if on_image_error is not None:
                on_image_error(i, total, pair, exc, tb_str)

    df_all_axons = pd.concat(all_axons, ignore_index=True) if all_axons else pd.DataFrame()
    df_all_images = pd.concat(all_images, ignore_index=True) if all_images else pd.DataFrame()
    df_all_mito = pd.concat(all_mito, ignore_index=True) if all_mito else pd.DataFrame()
    df_all_qc_report = pd.concat(all_qc_reports, ignore_index=True) if all_qc_reports else pd.DataFrame()

    return DatasetResult(
        df_axons=df_all_axons,
        df_image=df_all_images,
        df_mito=df_all_mito,
        df_qc_report=df_all_qc_report,
        errors=errors,
    )
