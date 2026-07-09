"""All disk reads of TEM/mask images go through this module."""

import numpy as np


def read_image(path) -> np.ndarray:
    """Read a TEM or mask TIFF/PNG/GIF, falling back to PIL if tifffile fails."""
    try:
        import tifffile as tiff
        arr = tiff.imread(path)
    except Exception:
        from PIL import Image
        arr = np.array(Image.open(path))
    if arr.ndim > 2:
        arr = arr[..., 0]
    return arr


def read_mask(path) -> np.ndarray:
    arr = read_image(path)
    if not np.issubdtype(arr.dtype, np.integer):
        arr = arr.astype(np.int64)
    return arr


def mask_value_histogram(mask: np.ndarray) -> dict:
    """Unique pixel values and counts. Powers the GUI's mask-info line."""
    unique, counts = np.unique(mask, return_counts=True)
    return {int(v): int(c) for v, c in zip(unique, counts)}
