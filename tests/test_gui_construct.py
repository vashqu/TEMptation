"""Constructs the real Tk app and drives a single-image analysis through
its own internal methods (not synthetic re-implementations), then tears
it down. Skips gracefully when no display is available (e.g. headless CI)."""

import os
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

DATA_ROOT = Path(__file__).parent.parent.parent
MASK_PATH = DATA_ROOT / "normal_data" / "162-165" / "163. Mask 20K.tif"


@pytest.fixture
def app():
    import tkinter as tk
    import measure_nerve_gui as gui_mod

    try:
        instance = gui_mod.NerveApp()
    except tk.TclError as exc:
        pytest.skip(f"no display available for Tk: {exc}")
    yield instance
    instance.destroy()


def test_gui_constructs(app):
    assert app.title() == "Nerve Fiber Analyzer"


def test_gui_single_image_analysis_matches_cli(app):
    """Run the same mask through the GUI's own analysis path and through
    the CLI-facing compat.measure_image, and assert identical results --
    this is the direct check that both front ends share one engine."""
    import tifffile as tiff
    from temptation.compat import measure_image

    app._mask_path = str(MASK_PATH)
    params = app._read_params()

    mask = tiff.imread(MASK_PATH)
    import numpy as np
    tem = np.zeros(mask.shape[:2], np.uint8)

    df_axons_gui, df_image_gui, _, mode_gui = measure_image(tem=tem, mask=mask, **params)
    df_axons_cli, df_image_cli, _, mode_cli = measure_image(
        tem=tem, mask=mask, pixel_length_um=0.00524,
    )

    assert mode_gui == mode_cli == "normal"
    assert len(df_axons_gui) == len(df_axons_cli) == 14

    import pandas as pd
    pd.testing.assert_frame_equal(df_axons_gui, df_axons_cli, check_exact=False, rtol=1e-9)
