#!/usr/bin/env python3

"""
measure_nerve.py
----------------
Segment and measure nerve fibers from TEM images + segmentation masks.

Usage examples
--------------
# Single image (normal, weighted watershed)
python measure_nerve.py --tem tem_001.tif --mask mask_001.tif \
    --bar 191 1 --plot

# Whole folder (auto-detect mode, simple watershed, custom labels)
python measure_nerve.py --folder ./data \
    --pixel-um 0.00524 --mode pathological --plot

# Force weighted watershed with area weight
python measure_nerve.py --tem img.tif --mask mask.tif \
    --bar 191 1 \
    --watershed-mode weighted --watershed-weight area \
    --plot
"""

import argparse
import re
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import ndimage as ndi
from scipy.stats import skew as scipy_skew
from skimage.measure import label, regionprops
from skimage.segmentation import watershed
from skimage.morphology import (
    binary_opening, binary_closing, binary_dilation,
    disk, remove_small_objects,
)
from scipy.spatial import cKDTree


# ---------------------------------------------------------------------------
# Helper: convexity (convex_perimeter / perimeter)
# ---------------------------------------------------------------------------

def _convexity(region) -> float:
    """
    Convexity = convex_hull_perimeter / actual_perimeter.
    Values close to 1 = convex shape; lower = irregular/concave boundary.
    Uses the regionprops convex_image to estimate convex perimeter.
    """
    from skimage.measure import perimeter as sk_perimeter
    convex_perim = sk_perimeter(region.convex_image)
    actual_perim = region.perimeter
    if actual_perim == 0:
        return np.nan
    return convex_perim / actual_perim


# ---------------------------------------------------------------------------
# Core segmentation & measurement
# ---------------------------------------------------------------------------

def _build_weighted_distance(axon_lab: np.ndarray, weight: str, beta: float) -> np.ndarray:
    """
    Return a weighted distance map for watershed.

    elevation = max(distance - beta * weight_map, 0)

    weight='radius' : weight_map = equivalent_radius  (larger axons expand more)
    weight='area'   : weight_map = sqrt(area)          (similar idea, area-based)
    beta            : bias strength; higher = larger axons claim more territory
    """
    distance, indices = ndi.distance_transform_edt(axon_lab == 0, return_indices=True)
    nearest_axon_map = axon_lab[indices[0], indices[1]]

    props = regionprops(axon_lab)
    max_label = axon_lab.max()
    weight_lookup = np.ones(max_label + 1, dtype=np.float32)   # default 1 (safe)

    for p in props:
        if weight == "radius":
            weight_lookup[p.label] = max(p.equivalent_diameter / 2.0, 1e-6)
        elif weight == "area":
            weight_lookup[p.label] = max(np.sqrt(p.area), 1e-6)
        else:
            raise ValueError(f"Unknown watershed weight: {weight!r}")

    weight_map = weight_lookup[nearest_axon_map]
    elevation = distance - beta * weight_map
    return np.maximum(elevation, 0.0)


def measure_image(
    tem: np.ndarray,
    mask: np.ndarray,
    pixel_length_um: float,
    # Label values
    myelin_val: int = 64,
    axoplasm_val: int = 192,
    mito_val: int = 128,
    # Mode
    mode: str = "auto",                        # "normal" | "pathological" | "auto"
    myelin_threshold_px: int = 200,
    # Morphological parameters
    smoothing_radius_px: int = 1,
    min_axon_area_px: int = 200,
    min_myelin_area_px: int = 300,
    # Watershed parameters
    watershed_mode: str = "weighted",          # "simple" | "weighted"
    watershed_weight: str = "radius",          # "radius" | "area"
    watershed_compactness: float = 0.001,
    watershed_beta: float = 1.0,
    # Detached myelin
    assign_detached_myelin: str = "none",      # "none" | "nearest"
):
    """
    Segment and measure nerve fibers.

    Returns
    -------
    df_axons      : pd.DataFrame  — one row per detected axon
    df_image      : pd.DataFrame  — one-row image-level summary
    labels_ws     : np.ndarray    — label image (matches axon_id in df_axons)
    resolved_mode : str           — "normal" or "pathological" (after auto-detect)

    New features added
    ------------------
    Axon-level:
      circularity        : 4π·area / perimeter²  (1 = perfect circle)
      convexity          : convex_hull_perimeter / perimeter
      axon_vol_fraction  : axon_area / fiber_area  (AVF)
      myelin_vol_fraction: myelin_area / fiber_area (MVF, NaN in pathological)

    Mitochondria (per axon, aggregated over individual mito objects):
      mito_mean_circularity   : mean of 4π·area/perimeter² per mito
      mito_mean_form_factor   : mean of perimeter²/(4π·area) per mito
      mito_mean_feret_um      : mean of max Feret diameter (µm) per mito
      mito_std_area_um2       : std of individual mito areas (size heterogeneity)
      mito_cv_area            : coefficient of variation of mito areas
      mito_area_skewness      : skewness of mito area distribution
      mito_max_area_um2       : largest single mito area in this axon
      mito_mean_dist_centroid_um : mean distance of mito centroids from axon centroid
      mito_std_dist_centroid_um  : std of those distances (spatial clustering)
    """
    px2 = pixel_length_um ** 2

    # ------------------------------------------------------------------
    # 1. Resolve mode (auto-detect if needed)
    # ------------------------------------------------------------------
    myelin_px_count = int((mask == myelin_val).sum())

    if mode == "auto":
        resolved_mode = "normal" if myelin_px_count >= myelin_threshold_px else "pathological"
    elif mode in ("normal", "pathological"):
        resolved_mode = mode
    else:
        raise ValueError(f"mode must be 'normal', 'pathological', or 'auto', got {mode!r}")

    # ------------------------------------------------------------------
    # 2. Binary masks & cleanup
    # ------------------------------------------------------------------
    axon_only = (mask == axoplasm_val)
    mito      = (mask == mito_val)
    myelin    = (mask == myelin_val)

    # Include mitochondria holes that fall inside axon territory
    axon_dil     = binary_dilation(axon_only, disk(1))
    mito_in_axon = mito & axon_dil
    axoplasm     = axon_only | mito_in_axon

    if smoothing_radius_px > 0:
        axoplasm = binary_opening(axoplasm, disk(smoothing_radius_px))
        myelin   = binary_closing(myelin,   disk(smoothing_radius_px))

    myelin = myelin & (~axoplasm)           # ensure mutual exclusion

    # ------------------------------------------------------------------
    # 3. Axon seeds
    # ------------------------------------------------------------------
    axon_lab = label(axoplasm, connectivity=2)
    axon_lab = remove_small_objects(axon_lab, min_size=min_axon_area_px)
    axon_lab = label(axon_lab > 0, connectivity=2)

    # ------------------------------------------------------------------
    # 4. Watershed
    # ------------------------------------------------------------------
    if resolved_mode == "normal":
        allowed_mask = (axon_lab > 0) | myelin
    else:
        allowed_mask = axon_lab > 0

    if watershed_mode == "weighted":
        elevation = _build_weighted_distance(axon_lab, watershed_weight, beta=watershed_beta)
    else:  # simple
        elevation = ndi.distance_transform_edt(axon_lab == 0)

    labels_ws = watershed(
        elevation,
        markers=axon_lab,
        mask=allowed_mask,
        compactness=watershed_compactness,
    )

    # ------------------------------------------------------------------
    # 4b. Nearest-axon map (pre-computed for detached myelin assignment)
    # ------------------------------------------------------------------
    if assign_detached_myelin == "nearest":
        _, _nn_idx = ndi.distance_transform_edt(axon_lab == 0, return_indices=True)
        nearest_axon_map_global = axon_lab[_nn_idx[0], _nn_idx[1]]
    else:
        nearest_axon_map_global = None

    # ------------------------------------------------------------------
    # 5. Metric extraction
    # ------------------------------------------------------------------
    H, W = mask.shape[:2]
    fov_area_um2 = H * W * px2

    props = regionprops(labels_ws)
    rows  = []

    for p in props:
        fiber_id = p.label
        min_row, min_col, max_row, max_col = p.bbox
        fiber_crop = p.image                          # boolean, bbox-sized

        axon_crop       = axoplasm[min_row:max_row, min_col:max_col]
        axon_mask_local = fiber_crop & axon_crop

        axon_area_px = int(axon_mask_local.sum())
        if axon_area_px == 0:
            continue

        if assign_detached_myelin == "nearest":
            myelin_area_px = int((myelin & (nearest_axon_map_global == fiber_id)).sum())
            fiber_area_px  = axon_area_px + myelin_area_px
        else:
            fiber_area_px  = p.area
            myelin_area_px = fiber_area_px - axon_area_px

        # --- Unit conversion ---
        fiber_area_um2  = fiber_area_px  * px2
        axon_area_um2   = axon_area_px   * px2
        myelin_area_um2 = myelin_area_px * px2

        # --- Geometric diameters ---
        d_inner = 2.0 * np.sqrt(axon_area_um2  / np.pi)
        d_outer = 2.0 * np.sqrt(fiber_area_um2 / np.pi)

        if resolved_mode == "normal" and d_outer > 0 and myelin_area_px >= min_myelin_area_px:
            g_ratio          = d_inner / d_outer
            myelin_thickness = (d_outer - d_inner) / 2.0
        else:
            g_ratio          = np.nan
            myelin_thickness = np.nan
            myelin_area_um2  = np.nan
            d_outer          = np.nan

        # --- Axon morphology (on the axon mask) ---
        axon_labeled_local = label(axon_mask_local)
        axon_rp = regionprops(axon_labeled_local)
        if not axon_rp:
            continue
        axon_props_local = max(axon_rp, key=lambda r: r.area)

        perimeter_um = axon_props_local.perimeter * pixel_length_um
        eccentricity = axon_props_local.eccentricity
        solidity     = axon_props_local.solidity

        # --- NEW: axon circularity ---
        # 4π·area / perimeter²  (pixel units, then converted)
        axon_perim_px = axon_props_local.perimeter
        if axon_perim_px > 0:
            circularity = (4.0 * np.pi * axon_area_px) / (axon_perim_px ** 2)
        else:
            circularity = np.nan

        # --- NEW: axon convexity (convex_hull_perimeter / perimeter) ---
        convexity = _convexity(axon_props_local)

        # --- NEW: axon volume fraction (AVF) and myelin volume fraction (MVF) ---
        avf = axon_area_um2 / fiber_area_um2 if fiber_area_um2 > 0 else np.nan
        if resolved_mode == "normal" and fiber_area_um2 > 0:
            mvf = myelin_area_um2 / fiber_area_um2
        else:
            mvf = np.nan

        cy_local, cx_local = axon_props_local.centroid
        cy_global = min_row + cy_local
        cx_global = min_col + cx_local
        cx_um = cx_global * pixel_length_um
        cy_um = cy_global * pixel_length_um

        # --- Mitochondria ---
        mito_crop         = mito[min_row:max_row, min_col:max_col]
        mito_in_axon_crop = mito_crop & axon_mask_local
        mito_area_px_loc  = int(mito_in_axon_crop.sum())
        mito_area_um2     = mito_area_px_loc * px2
        mito_lab          = label(mito_in_axon_crop, connectivity=2)
        mito_regions      = regionprops(mito_lab)
        mito_count        = len(mito_regions)
        mito_density      = mito_count / axon_area_um2 if axon_area_um2 > 0 else np.nan

        # --- NEW: per-mito shape features, aggregated per axon ---
        if mito_count > 0:
            mito_circularities = []
            mito_form_factors  = []
            mito_ferets_um     = []
            mito_areas_um2     = []
            mito_cx_list       = []   # centroid x in local crop coords
            mito_cy_list       = []

            for mr in mito_regions:
                m_area = mr.area
                m_perim = mr.perimeter
                m_areas_um2 = m_area * px2
                mito_areas_um2.append(m_areas_um2)

                # circularity
                if m_perim > 0:
                    m_circ = (4.0 * np.pi * m_area) / (m_perim ** 2)
                else:
                    m_circ = np.nan
                mito_circularities.append(m_circ)

                # form factor (inverse of circularity — > 1 means elongated/complex)
                if m_perim > 0 and m_area > 0:
                    m_ff = (m_perim ** 2) / (4.0 * np.pi * m_area)
                else:
                    m_ff = np.nan
                mito_form_factors.append(m_ff)

                # Feret's diameter (max caliper) — in pixels, convert to µm
                m_feret_um = mr.feret_diameter_max * pixel_length_um
                mito_ferets_um.append(m_feret_um)

                # centroid in local crop coords
                mito_cy_list.append(mr.centroid[0])
                mito_cx_list.append(mr.centroid[1])

            mito_areas_arr = np.array(mito_areas_um2)

            mito_mean_circularity = float(np.nanmean(mito_circularities))
            mito_mean_form_factor = float(np.nanmean(mito_form_factors))
            mito_mean_feret_um    = float(np.nanmean(mito_ferets_um))
            mito_std_area_um2     = float(np.std(mito_areas_arr)) if mito_count > 1 else 0.0
            mito_max_area_um2     = float(np.max(mito_areas_arr))
            mito_cv_area          = (
                mito_std_area_um2 / float(np.mean(mito_areas_arr))
                if mito_count > 1 and np.mean(mito_areas_arr) > 0
                else np.nan
            )
            mito_area_skewness    = (
                float(scipy_skew(mito_areas_arr))
                if mito_count > 2
                else np.nan
            )

            # spatial: distances of mito centroids from axon centroid (local coords)
            mito_cx_arr = np.array(mito_cx_list)
            mito_cy_arr = np.array(mito_cy_list)
            dists_px = np.sqrt((mito_cx_arr - cx_local) ** 2 +
                               (mito_cy_arr - cy_local) ** 2)
            dists_um = dists_px * pixel_length_um
            mito_mean_dist_centroid_um = float(np.mean(dists_um))
            mito_std_dist_centroid_um  = float(np.std(dists_um)) if mito_count > 1 else 0.0

        else:
            mito_mean_circularity      = np.nan
            mito_mean_form_factor      = np.nan
            mito_mean_feret_um         = np.nan
            mito_std_area_um2          = np.nan
            mito_max_area_um2          = np.nan
            mito_cv_area               = np.nan
            mito_area_skewness         = np.nan
            mito_mean_dist_centroid_um = np.nan
            mito_std_dist_centroid_um  = np.nan

        rows.append({
            # --- identifiers ---
            "axon_id":                     fiber_id,
            # --- areas ---
            "axon_area_um2":               axon_area_um2,
            "myelin_area_um2":             myelin_area_um2,   # NaN in pathological
            "fiber_area_um2":              fiber_area_um2,
            # --- diameters & myelin ---
            "d_inner_um":                  d_inner,
            "d_outer_um":                  d_outer,           # NaN in pathological
            "g_ratio":                     g_ratio,           # NaN in pathological
            "myelin_thickness_um":         myelin_thickness,  # NaN in pathological
            # --- axon shape ---
            "perimeter_um":                perimeter_um,
            "eccentricity":                eccentricity,
            "solidity":                    solidity,
            "circularity":                 circularity,       # NEW
            "convexity":                   convexity,         # NEW
            # --- volume fractions ---
            "axon_vol_fraction":           avf,               # NEW (AVF)
            "myelin_vol_fraction":         mvf,               # NEW (MVF, NaN in pathological)
            # --- centroids ---
            "centroid_x_um":               cx_um,
            "centroid_y_um":               cy_um,
            "centroid_x_px":               cx_global,
            "centroid_y_px":               cy_global,
            # --- mitochondria counts / bulk ---
            "mito_count":                  mito_count,
            "mito_area_um2":               mito_area_um2,
            "mito_density_per_um2":        mito_density,
            # --- mitochondria shape (NEW) ---
            "mito_mean_circularity":       mito_mean_circularity,
            "mito_mean_form_factor":       mito_mean_form_factor,
            "mito_mean_feret_um":          mito_mean_feret_um,
            "mito_std_area_um2":           mito_std_area_um2,
            "mito_max_area_um2":           mito_max_area_um2,
            "mito_cv_area":                mito_cv_area,
            "mito_area_skewness":          mito_area_skewness,
            # --- mitochondria spatial (NEW) ---
            "mito_mean_dist_centroid_um":  mito_mean_dist_centroid_um,
            "mito_std_dist_centroid_um":   mito_std_dist_centroid_um,
        })

    df_axons = pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # 6. Nearest-neighbour distances
    # ------------------------------------------------------------------
    if not df_axons.empty:
        if len(df_axons) >= 2:
            pts = df_axons[["centroid_x_um", "centroid_y_um"]].to_numpy()
            tree = cKDTree(pts)
            dists, _ = tree.query(pts, k=2)
            df_axons["nearest_neighbor_um"] = dists[:, 1]
        else:
            df_axons["nearest_neighbor_um"] = np.nan

    # ------------------------------------------------------------------
    # 7. Image-level summary
    # ------------------------------------------------------------------
    if not df_axons.empty:
        n_fibers = len(df_axons)
        fiber_density_per_mm2 = (n_fibers / fov_area_um2) * 1e6

        summary = {
            "n_fibers":              n_fibers,
            "fov_area_um2":          fov_area_um2,
            "fiber_density_per_mm2": fiber_density_per_mm2,
            "mean_axon_area_um2":    df_axons["axon_area_um2"].mean(),
            "mean_circularity":      df_axons["circularity"].mean(),
            "mean_convexity":        df_axons["convexity"].mean(),
            "mean_avf":              df_axons["axon_vol_fraction"].mean(),
            "mode":                  resolved_mode,
        }
        if resolved_mode == "normal":
            summary["mean_g_ratio"]             = df_axons["g_ratio"].mean()
            summary["mean_myelin_thickness_um"] = df_axons["myelin_thickness_um"].mean()
            summary["mean_mvf"]                 = df_axons["myelin_vol_fraction"].mean()
        else:
            summary["mean_g_ratio"]             = np.nan
            summary["mean_myelin_thickness_um"] = np.nan
            summary["mean_mvf"]                 = np.nan

        mito_outside          = (mask == mito_val) & ~(axon_only | myelin)
        mito_outside_area_um2 = float(mito_outside.sum()) * px2
        total_mito_px         = int((mask == mito_val).sum())
        mito_outside_frac     = float(mito_outside.sum()) / max(total_mito_px, 1)
        summary["mito_outside_area_um2"] = mito_outside_area_um2
        summary["mito_outside_frac"]     = mito_outside_frac

        total_myelin_px             = int((mask == myelin_val).sum())
        total_myelin_um2            = total_myelin_px * px2
        myelin_area_fraction_of_fov = total_myelin_um2 / fov_area_um2 if fov_area_um2 > 0 else 0.0
        myelin_lab                  = label(mask == myelin_val)
        myelin_component_count      = int(myelin_lab.max())
        summary["total_myelin_area_um2"]       = total_myelin_um2
        summary["myelin_area_fraction_of_fov"] = myelin_area_fraction_of_fov
        summary["myelin_component_count"]      = myelin_component_count

        df_image = pd.DataFrame([summary])
    else:
        df_image = pd.DataFrame()

    return df_axons, df_image, labels_ws, resolved_mode


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plot(
    tem: np.ndarray,
    mask: np.ndarray,
    labels_ws: np.ndarray,
    df_axons: pd.DataFrame,
    myelin_val: int = 64,
    axoplasm_val: int = 192,
    mito_val: int = 128,
    title: str = "",
    save_path: Path = None,
):
    """
    2-panel figure:
      Left  — raw TEM
      Right — TEM + tissue-type overlay (myelin / axoplasm / mito)
               with axon ID numbers at each axon centroid.

    Overlay colours (fixed, tissue-type based):
      Myelin    — steel blue,  alpha 0.45
      Axoplasm  — orange,      alpha 0.50
      Mito      — red,         alpha 0.80
    """
    MYELIN_RGBA   = np.array([0.27, 0.51, 0.71, 0.45], dtype=np.float32)
    AXOPLASM_RGBA = np.array([1.00, 0.60, 0.10, 0.50], dtype=np.float32)
    MITO_RGBA     = np.array([0.85, 0.15, 0.15, 0.80], dtype=np.float32)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(title, fontsize=13)

    axes[0].imshow(tem, cmap="gray", interpolation="nearest")
    axes[0].set_title("TEM image")
    axes[0].axis("off")

    axes[1].imshow(tem, cmap="gray", interpolation="nearest")

    overlay = np.zeros((*labels_ws.shape, 4), dtype=np.float32)

    # Tissue overlays are based solely on raw mask values, NOT restricted by
    # watershed labels_ws.  labels_ws is used only for contours and axon IDs.
    axon_only     = (mask == axoplasm_val)
    mito_raw      = (mask == mito_val)
    axon_dil      = binary_dilation(axon_only, disk(1))
    mito_in_axon  = mito_raw & axon_dil
    axoplasm_full = axon_only | mito_in_axon

    myelin_mask   = (mask == myelin_val)
    axoplasm_mask = axoplasm_full.copy()
    mito_mask     = ndi.binary_fill_holes(mito_raw)

    axoplasm_mask = axoplasm_mask & ~mito_mask      # red sits on top of orange

    overlay[myelin_mask]   = MYELIN_RGBA
    overlay[axoplasm_mask] = AXOPLASM_RGBA
    overlay[mito_mask]     = MITO_RGBA

    axes[1].imshow(overlay, interpolation="nearest")

    from skimage import measure as sk_measure
    unique_labels = np.unique(labels_ws)
    unique_labels = unique_labels[unique_labels > 0]
    for lbl in unique_labels:
        fiber_bin = (labels_ws == lbl).astype(np.uint8)
        contours = sk_measure.find_contours(fiber_bin, level=0.5)
        for contour in contours:
            axes[1].plot(
                contour[:, 1], contour[:, 0],
                color="white", linewidth=0.8, alpha=0.9,
            )

    if not df_axons.empty:
        for _, row in df_axons.iterrows():
            axes[1].text(
                row["centroid_x_px"], row["centroid_y_px"],
                str(int(row["axon_id"])),
                color="white", fontsize=6, ha="center", va="center",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.45, lw=0),
            )

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor=MYELIN_RGBA[:3],   alpha=0.8, label="Myelin"),
        Patch(facecolor=AXOPLASM_RGBA[:3], alpha=0.8, label="Axoplasm"),
        Patch(facecolor=MITO_RGBA[:3],     alpha=0.8, label="Mitochondria"),
        Line2D([0], [0], color="white", linewidth=1.2, label="Fiber boundary"),
    ]
    axes[1].legend(handles=legend_elements, loc="lower right",
                   fontsize=7, framealpha=0.6)

    axes[1].set_title("Overlay + axon IDs")
    axes[1].axis("off")

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"  [plot] saved → {save_path}")

    plt.show(block=False)
    plt.pause(0.1)
    plt.close(fig)


# ---------------------------------------------------------------------------
# File discovery helpers
# ---------------------------------------------------------------------------

def _find_pairs_in_folder(folder: Path):
    """
    Scan folder for paired TEM + mask files.

    Supported naming conventions
    ----------------------------
    Pattern A — explicit prefix with underscore separator:
        tem_<id>.tif   OR   axon_<id>.tif   →  TEM image
        mask_<id>.tif                        →  mask image

    Pattern B — numeric prefix with tissue keyword anywhere in name:
        <number>. Axon *.tif  OR  <number>. TEM *.tif  →  TEM image
        <number>. Mask *.tif                            →  mask image
        IDs are matched by the leading integer.

    Returns list of dicts: {id, tem_path, mask_path}
    """
    tif_files = sorted(folder.glob("*.tif")) + sorted(folder.glob("*.tiff"))

    tem_map  = {}
    mask_map = {}

    for f in tif_files:
        stem       = f.stem
        stem_lower = stem.lower()

        # Pattern A: tem_<id> / axon_<id> / mask_<id>
        if stem_lower.startswith("tem_") or stem_lower.startswith("axon_"):
            img_id = stem.split("_", 1)[1]
            tem_map.setdefault(img_id, f)
        elif stem_lower.startswith("mask_"):
            img_id = stem.split("_", 1)[1]
            mask_map.setdefault(img_id, f)

        # Pattern B: leading number + tissue keyword anywhere in name
        else:
            m = re.match(r"^(\d+)", stem)
            if not m:
                continue
            img_id = m.group(1)
            if "mask" in stem_lower:
                mask_map.setdefault(img_id, f)
            elif "axon" in stem_lower or "tem" in stem_lower:
                tem_map.setdefault(img_id, f)

    def _sort_key(x):
        return int(x) if x.isdigit() else x

    pairs = []
    all_ids = set(tem_map) | set(mask_map)
    for img_id in sorted(all_ids, key=_sort_key):
        if img_id not in tem_map:
            warnings.warn(f"[SKIP] No TEM file for id={img_id!r}")
            continue
        if img_id not in mask_map:
            warnings.warn(f"[SKIP] No mask file for id={img_id!r}")
            continue
        pairs.append({
            "id":        img_id,
            "tem_path":  tem_map[img_id],
            "mask_path": mask_map[img_id],
        })

    return pairs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Segment and measure nerve fibers from TEM + mask images.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    input_grp = p.add_argument_group("Input (use --folder OR --tem + --mask)")
    input_grp.add_argument("--folder", type=Path, default=None)
    input_grp.add_argument("--tem",    type=Path, default=None)
    input_grp.add_argument("--mask",   type=Path, default=None)

    scale_grp = p.add_argument_group("Scale (provide exactly one option)")
    scale_ex = scale_grp.add_mutually_exclusive_group()
    scale_ex.add_argument("--pixel-um", type=float, default=None)
    scale_ex.add_argument("--bar", nargs=2, metavar=("PX", "UM"), type=float, default=None)

    p.add_argument("--mode", choices=["normal", "pathological", "auto"], default="auto")
    p.add_argument("--myelin-threshold", type=int, default=200)
    p.add_argument("--myelin-val",   type=int, default=64)
    p.add_argument("--axoplasm-val", type=int, default=192)
    p.add_argument("--mito-val",     type=int, default=128)
    p.add_argument("--smoothing-radius", type=int,   default=1)
    p.add_argument("--min-axon-area",    type=int,   default=200)
    p.add_argument("--min-myelin-area",  type=int,   default=300)
    p.add_argument("--watershed-mode",   choices=["simple", "weighted"], default="weighted")
    p.add_argument("--watershed-weight", choices=["radius", "area"],     default="radius")
    p.add_argument("--watershed-compactness", type=float, default=0.001)
    p.add_argument("--watershed-beta",        type=float, default=1.0,
                   help="Bias strength for weighted watershed; higher = larger axons claim more territory.")
    p.add_argument("--assign-detached-myelin", choices=["none", "nearest"], default="none",
                   help="Assign detached myelin pixels to the nearest axon for myelin metric computation.")
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--plot", action="store_true")

    return p


def resolve_pixel_length(args) -> float:
    if args.pixel_um is not None:
        return args.pixel_um
    elif args.bar is not None:
        bar_px, bar_um = args.bar
        return bar_um / bar_px
    else:
        raise SystemExit(
            "ERROR: You must provide either --pixel-um or --bar <PX> <UM> to set the scale."
        )


def process_pair(
    image_id: str,
    tem_path: Path,
    mask_path: Path,
    pixel_length_um: float,
    args,
    output_dir: Path,
) -> tuple:
    tem  = tiff.imread(tem_path)
    mask = tiff.imread(mask_path)

    df_axons, df_image, labels_ws, resolved_mode = measure_image(
        tem=tem,
        mask=mask,
        pixel_length_um=pixel_length_um,
        myelin_val=args.myelin_val,
        axoplasm_val=args.axoplasm_val,
        mito_val=args.mito_val,
        mode=args.mode,
        myelin_threshold_px=args.myelin_threshold,
        smoothing_radius_px=args.smoothing_radius,
        min_axon_area_px=args.min_axon_area,
        min_myelin_area_px=args.min_myelin_area,
        watershed_mode=args.watershed_mode,
        watershed_weight=args.watershed_weight,
        watershed_compactness=args.watershed_compactness,
        watershed_beta=args.watershed_beta,
        assign_detached_myelin=args.assign_detached_myelin,
    )

    print(f"  [{resolved_mode.upper()}] id={image_id!r}  "
          f"→ {len(df_axons)} axons detected")

    for df in (df_axons, df_image):
        if not df.empty:
            if "image_id" not in df.columns:
                df.insert(0, "image_id", image_id)
            else:
                df["image_id"] = image_id
            if "mode" not in df.columns:
                df.insert(1, "mode", resolved_mode)
            else:
                df["mode"] = resolved_mode

    if args.plot and not df_axons.empty:
        plot_path = output_dir / f"overlay_{image_id}.png"
        try:
            make_plot(
                tem=tem,
                mask=mask,
                labels_ws=labels_ws,
                df_axons=df_axons,
                myelin_val=args.myelin_val,
                axoplasm_val=args.axoplasm_val,
                mito_val=args.mito_val,
                title=f"Image {image_id}  [{resolved_mode}]",
                save_path=plot_path,
            )
        except Exception as plot_exc:
            print(f"  [WARNING] Plot failed for id={image_id!r}: {plot_exc}")
            traceback.print_exc()

    return df_axons, df_image


def main():
    parser = build_parser()
    args   = parser.parse_args()

    if args.folder is not None:
        if args.tem is not None or args.mask is not None:
            parser.error("Use either --folder OR (--tem + --mask), not both.")
        if not args.folder.is_dir():
            parser.error(f"--folder {args.folder!r} is not a directory.")
        pairs = _find_pairs_in_folder(args.folder)
        if not pairs:
            sys.exit(f"ERROR: No valid TEM+mask pairs found in {args.folder}")
        output_dir = args.output_dir or args.folder
    else:
        if args.tem is None or args.mask is None:
            parser.error("Provide --folder OR both --tem and --mask.")
        if not args.tem.is_file():
            parser.error(f"--tem file not found: {args.tem}")
        if not args.mask.is_file():
            parser.error(f"--mask file not found: {args.mask}")
        pairs = [{"id": args.tem.stem, "tem_path": args.tem, "mask_path": args.mask}]
        output_dir = args.output_dir or args.tem.parent

    output_dir.mkdir(parents=True, exist_ok=True)
    pixel_length_um = resolve_pixel_length(args)
    print(f"Pixel length: {pixel_length_um:.6f} µm/px")

    all_axons  = []
    all_images = []

    for pair in pairs:
        try:
            df_axons, df_image = process_pair(
                image_id=pair["id"],
                tem_path=pair["tem_path"],
                mask_path=pair["mask_path"],
                pixel_length_um=pixel_length_um,
                args=args,
                output_dir=output_dir,
            )
            if not df_axons.empty:
                all_axons.append(df_axons)
            if not df_image.empty:
                all_images.append(df_image)
        except Exception as exc:
            print(f"\n[ERROR] id={pair['id']!r} failed with: {exc}")
            traceback.print_exc()
            print()

    if all_axons:
        df_all_axons = pd.concat(all_axons, ignore_index=True)
        axon_csv = output_dir / "axons.csv"
        df_all_axons.to_csv(axon_csv, index=False)
        print(f"\nAxon-level results  → {axon_csv}  ({len(df_all_axons)} rows)")
    else:
        print("\nNo axons detected across all images.")
        df_all_axons = pd.DataFrame()

    if all_images:
        df_all_images = pd.concat(all_images, ignore_index=True)
        img_csv = output_dir / "image_summary.csv"
        df_all_images.to_csv(img_csv, index=False)
        print(f"Image-level summary → {img_csv}  ({len(df_all_images)} rows)")

    return df_all_axons


if __name__ == "__main__":
    main()