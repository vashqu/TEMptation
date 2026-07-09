"""Per-axon geometry and shape metrics."""

import numpy as np


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


def axon_geometry(
    axon_area_px: int,
    fiber_area_px: int,
    myelin_area_px: int,
    resolved_mode: str,
    min_myelin_area_px: int,
    pixel_length_um: float,
) -> dict:
    """Areas, diameters, g-ratio, myelin thickness, AVF/MVF.

    Verbatim port of measure_nerve.py:249-266, 289-294.
    """
    px2 = pixel_length_um ** 2

    fiber_area_um2 = fiber_area_px * px2
    axon_area_um2 = axon_area_px * px2
    myelin_area_um2 = myelin_area_px * px2

    d_inner = 2.0 * np.sqrt(axon_area_um2 / np.pi)
    d_outer = 2.0 * np.sqrt(fiber_area_um2 / np.pi)

    if resolved_mode == "normal" and d_outer > 0 and myelin_area_px >= min_myelin_area_px:
        g_ratio = d_inner / d_outer
        myelin_thickness = (d_outer - d_inner) / 2.0
    else:
        g_ratio = np.nan
        myelin_thickness = np.nan
        myelin_area_um2 = np.nan
        d_outer = np.nan

    avf = axon_area_um2 / fiber_area_um2 if fiber_area_um2 > 0 else np.nan
    if resolved_mode == "normal" and fiber_area_um2 > 0:
        mvf = myelin_area_um2 / fiber_area_um2
    else:
        mvf = np.nan

    return {
        "axon_area_um2": axon_area_um2,
        "myelin_area_um2": myelin_area_um2,
        "fiber_area_um2": fiber_area_um2,
        "d_inner_um": d_inner,
        "d_outer_um": d_outer,
        "g_ratio": g_ratio,
        "myelin_thickness_um": myelin_thickness,
        "axon_vol_fraction": avf,
        "myelin_vol_fraction": mvf,
    }


def shape_metrics(axon_props_local, axon_area_px: int, pixel_length_um: float) -> dict:
    """perimeter, eccentricity, solidity, circularity, convexity.

    Verbatim port of measure_nerve.py:274-287. `axon_area_px` must be the
    full local-crop axon area (axon_mask_local.sum()), not
    axon_props_local.area -- these can differ when the crop has more than
    one connected fragment, and the original code uses the crop sum.
    """
    perimeter_um = axon_props_local.perimeter * pixel_length_um
    eccentricity = axon_props_local.eccentricity
    solidity = axon_props_local.solidity

    axon_perim_px = axon_props_local.perimeter
    if axon_perim_px > 0:
        circularity = (4.0 * np.pi * axon_area_px) / (axon_perim_px ** 2)
    else:
        circularity = np.nan

    convexity = _convexity(axon_props_local)

    return {
        "perimeter_um": perimeter_um,
        "eccentricity": eccentricity,
        "solidity": solidity,
        "circularity": circularity,
        "convexity": convexity,
    }
