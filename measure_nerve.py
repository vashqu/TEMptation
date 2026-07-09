#!/usr/bin/env python3

"""
measure_nerve.py
----------------
Segment and measure nerve fibers from TEM images + segmentation masks.

This file is now a thin shim over the `temptation` package (see
IMPLEMENTATION_BLUEPRINT.md). All engine logic lives in temptation/*.py;
this module exists so existing imports and CLI invocations keep working
unchanged:

    from measure_nerve import measure_image, _find_pairs_in_folder

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

import sys
from pathlib import Path

_this_dir = Path(__file__).resolve().parent
if str(_this_dir) not in sys.path:
    sys.path.insert(0, str(_this_dir))

from temptation.compat import measure_image  # noqa: F401,E402
from temptation.discovery import find_pairs as _find_pairs_in_folder  # noqa: F401,E402
from temptation.plotting import overlay_figure  # noqa: F401,E402
from temptation.cli import (  # noqa: F401,E402
    build_parser,
    resolve_pixel_length,
    process_pair,
    main,
)


def make_plot(tem, mask, labels_ws, df_axons, myelin_val=64, axoplasm_val=192,
              mito_val=128, title="", save_path=None):
    """Back-compat wrapper over temptation.plotting; matches the original
    make_plot() signature and behavior (save + non-blocking show)."""
    from temptation import plotting
    fig = plotting.overlay_figure(
        tem, mask, labels_ws, df_axons,
        myelin_val=myelin_val, axoplasm_val=axoplasm_val, mito_val=mito_val,
        title=title,
    )
    if save_path is not None:
        plotting.save_overlay(fig, save_path)
    plotting.show_overlay_nonblocking(fig)


if __name__ == "__main__":
    main()
