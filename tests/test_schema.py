"""Schema stability tests (Phase 6, fixes F8's "unstable schema when 0
rows" half). IMPLEMENTATION_BLUEPRINT.md Sec 9 Phase 6 Checks: an
all-background mask yields axons.csv with 0 rows and the full header; a
1-axon image emits nearest_neighbor_um = NaN rather than omitting the
column."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from temptation.cli import main as cli_main
from temptation.config import SegmentationConfig
from temptation.pipeline import analyze_image_legacy
from temptation.schema import AXON_COLUMNS_V2, IDENTITY_COLUMNS, IMAGE_COLUMNS_V2, MITO_COLUMNS

# cli.py appends image_id/mode (already in AXON_COLUMNS_V2/IMAGE_COLUMNS_V2)
# plus IDENTITY_COLUMNS (group/image_path/mask_path/pixel_size_um/
# schema_version, not part of either V2 tuple) to every written row.
_WRITTEN_AXON_COLUMN_COUNT = len(AXON_COLUMNS_V2) + len(IDENTITY_COLUMNS)
_WRITTEN_IMAGE_COLUMN_COUNT = len(IMAGE_COLUMNS_V2) + len(IDENTITY_COLUMNS)

DATA_ROOT = Path(__file__).parent.parent.parent


def test_all_background_mask_zero_rows_full_header():
    mask = np.zeros((200, 200), dtype=np.uint8)
    tem = np.zeros((200, 200), dtype=np.uint8)
    cfg = SegmentationConfig()

    df_axons, df_image, labels_ws, mode, df_mito = analyze_image_legacy(tem, mask, 0.00524, cfg)

    assert df_axons.empty
    assert len(df_axons.columns) == len(AXON_COLUMNS_V2) - 2  # minus image_id, mode (added by cli.py)
    assert df_image.empty
    assert len(df_image.columns) == len(IMAGE_COLUMNS_V2) - 1  # minus image_id
    assert df_mito.empty
    assert len(df_mito.columns) == len(MITO_COLUMNS)


def test_single_axon_image_has_nearest_neighbor_um_column_as_nan():
    """A single-fiber image can't have a nearest-neighbor distance, but
    the column must still be present (NaN), not omitted."""
    from skimage.draw import disk as skdisk

    mask = np.zeros((300, 300), dtype=np.uint8)
    rr, cc = skdisk((150, 150), 60, shape=mask.shape)
    mask[rr, cc] = 64  # myelin ring
    rr, cc = skdisk((150, 150), 50, shape=mask.shape)
    mask[rr, cc] = 192  # axoplasm
    tem = np.zeros(mask.shape, dtype=np.uint8)
    cfg = SegmentationConfig()

    df_axons, df_image, labels_ws, mode, df_mito = analyze_image_legacy(tem, mask, 0.00524, cfg)

    assert len(df_axons) == 1
    assert "nearest_neighbor_um" in df_axons.columns
    assert np.isnan(df_axons["nearest_neighbor_um"].iloc[0])


def test_cli_writes_full_header_when_batch_entirely_empty(tmp_path, monkeypatch):
    """An all-background-mask 'dataset' (one image, zero axons) must
    still produce axons.csv/image_summary.csv with the correct header,
    not an empty (0-column) file."""
    import tifffile as tiff

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    mask = np.zeros((100, 100), dtype=np.uint8)
    tem = np.zeros((100, 100), dtype=np.uint8)
    tiff.imwrite(data_dir / "1. Axon 20K.tif", tem)
    tiff.imwrite(data_dir / "1. Mask 20K.tif", mask)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(data_dir), "--pixel-um", "0.00524",
        "--output-dir", str(out_dir),
    ])
    cli_main()

    axons = pd.read_csv(out_dir / "axons.csv")
    images = pd.read_csv(out_dir / "image_summary.csv")
    assert len(axons) == 0
    assert len(axons.columns) == _WRITTEN_AXON_COLUMN_COUNT
    assert len(images) == 0
    assert len(images.columns) == _WRITTEN_IMAGE_COLUMN_COUNT


def test_mixed_batch_empty_and_nonempty_images_consistent_schema(tmp_path, monkeypatch):
    """A batch with one empty image and one real image must still yield
    a single, consistent column set across both."""
    import shutil
    import tifffile as tiff

    data_dir = tmp_path / "data"
    data_dir.mkdir()

    # empty image
    tiff.imwrite(data_dir / "1. Axon 20K.tif", np.zeros((100, 100), dtype=np.uint8))
    tiff.imwrite(data_dir / "1. Mask 20K.tif", np.zeros((100, 100), dtype=np.uint8))

    # real image, copied from the actual dataset
    real_axon = DATA_ROOT / "normal_data" / "162-165" / "162. Axon 20K.tif"
    real_mask = DATA_ROOT / "normal_data" / "162-165" / "162. Mask 20K.tif"
    shutil.copy(real_axon, data_dir / "2. Axon 20K.tif")
    shutil.copy(real_mask, data_dir / "2. Mask 20K.tif")

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.setattr("sys.argv", [
        "prog", "--folder", str(data_dir), "--pixel-um", "0.00524",
        "--output-dir", str(out_dir),
    ])
    cli_main()

    axons = pd.read_csv(out_dir / "axons.csv")
    assert len(axons.columns) == _WRITTEN_AXON_COLUMN_COUNT
    assert set(axons["image_id"].astype(str)) == {"2"}  # image 1 contributed 0 rows, correctly
    assert len(axons) > 0
