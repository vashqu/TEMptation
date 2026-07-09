"""All CSV writes go through this module -- called identically by CLI and
GUI so file naming never drifts between the two front-ends."""

from pathlib import Path

import pandas as pd


def write_axon_csv(df_axons: pd.DataFrame, out_dir: Path) -> Path:
    path = Path(out_dir) / "axons.csv"
    df_axons.to_csv(path, index=False)
    return path


def write_image_csv(df_images: pd.DataFrame, out_dir: Path) -> Path:
    path = Path(out_dir) / "image_summary.csv"
    df_images.to_csv(path, index=False)
    return path
