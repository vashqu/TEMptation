"""Safe-math primitives shared by every metrics module.

Centralizing these means the NaN policy (CLAUDE.md Sec 5.1 / 12.2) is
defined exactly once instead of re-implemented ad hoc at each call site.
"""

import numpy as np


def safe_divide(numerator, denominator):
    """Return numerator/denominator, or NaN if denominator is 0/None/NaN."""
    if denominator is None:
        return np.nan
    try:
        if denominator == 0 or np.isnan(denominator):
            return np.nan
    except TypeError:
        return np.nan
    return numerator / denominator


def safe_cv(values, ddof=1):
    """Coefficient of variation (std/mean). NaN on empty input or zero mean."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.nan
    mean = np.nanmean(values)
    if mean == 0 or np.isnan(mean):
        return np.nan
    with np.errstate(invalid="ignore"):
        std = np.nanstd(values, ddof=ddof)
    return std / mean


def distribution_stats(values, prefix, ddof=1):
    """mean/std/cv/median/min/max/iqr over `values`, keyed as f"{prefix}_<stat>".

    Empty input -> every stat is NaN. Never raises.
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if values.size == 0:
        return {
            f"{prefix}_mean": np.nan,
            f"{prefix}_std": np.nan,
            f"{prefix}_cv": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_min": np.nan,
            f"{prefix}_max": np.nan,
            f"{prefix}_iqr": np.nan,
        }

    mean = float(np.mean(values))
    if values.size > ddof:
        with np.errstate(invalid="ignore"):
            std = float(np.std(values, ddof=ddof))
    else:
        std = np.nan
    cv = std / mean if (mean != 0 and not np.isnan(std)) else np.nan
    q75, q25 = np.percentile(values, [75, 25])

    return {
        f"{prefix}_mean": mean,
        f"{prefix}_std": std,
        f"{prefix}_cv": cv,
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_min": float(np.min(values)),
        f"{prefix}_max": float(np.max(values)),
        f"{prefix}_iqr": float(q75 - q25),
    }


def safe_aspect_ratio(region):
    """major_axis_length / minor_axis_length. NaN if minor_axis_length <= 0."""
    minor = getattr(region, "minor_axis_length", None)
    if minor is None or minor <= 0:
        return np.nan
    return region.major_axis_length / minor


def equivalent_radius(area):
    """sqrt(area / pi). NaN on negative/NaN/None input."""
    if area is None:
        return np.nan
    if np.isnan(area) or area < 0:
        return np.nan
    return float(np.sqrt(area / np.pi))
