#!/usr/bin/env python3

"""
AxonDeepSeg bridge for TEMptation.

This module is intentionally optional: importing it does not import AxonDeepSeg,
torch, nnU-Net, or napari. The heavy dependency stack is only touched inside a
subprocess when automatic segmentation is requested.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
import tifffile as tiff


ADS_MODEL_DESCRIPTIONS = {
    "generalist": (
        "Multi-domain axon and myelin segmentation model trained on TEM, SEM, "
        "BF and CARS data."
    ),
    "unmyelinated-TEM": (
        "Segments both myelinated and unmyelinated axons on TEM data. Also "
        "segments oligodendrocytes (nuclei and cytoplasmic processes)."
    ),
}

ADS_MODEL_CHOICES = tuple(ADS_MODEL_DESCRIPTIONS.keys())

MITOCHONDRIA_NOTE = (
    "Automatic AxonDeepSeg masks do not include mitochondria at this point; "
    "mitochondria metrics require manual mask refinement with the mitochondria "
    "label."
)


class AxonDeepSegError(RuntimeError):
    """Raised when AxonDeepSeg cannot be run or does not produce usable masks."""


@dataclass(frozen=True)
class AutoSegmentationResult:
    image_path: Path
    mask_path: Path
    model_name: str
    class_mask_paths: dict[str, Path]
    output_dir: Path
    ads_work_dir: Path | None
    mitochondria_note: str = MITOCHONDRIA_NOTE


def _default_ads_package_dir() -> Path | None:
    env_path = os.environ.get("AXONDEEPSEG_PATH")
    if env_path:
        return Path(env_path).expanduser()

    workspace_root = Path(__file__).resolve().parent.parent
    for dirname in ("axondeepseg", "AxonDeepSeg"):
        workspace_sibling = workspace_root / dirname
        if workspace_sibling.exists():
            return workspace_sibling

    return None


def _resolve_ads_python(ads_python: str | Path | None) -> str:
    if ads_python:
        return str(Path(ads_python).expanduser())
    env_python = os.environ.get("AXONDEEPSEG_PYTHON")
    if env_python:
        return env_python
    return sys.executable


def _build_ads_env(ads_package_dir: str | Path | None) -> dict[str, str]:
    env = os.environ.copy()
    package_dir = Path(ads_package_dir).expanduser() if ads_package_dir else _default_ads_package_dir()
    if package_dir and package_dir.exists():
        existing = env.get("PYTHONPATH")
        env["PYTHONPATH"] = (
            str(package_dir) if not existing else str(package_dir) + os.pathsep + existing
        )
    return env


def _tail(text: str, max_chars: int = 5000) -> str:
    text = text.strip()
    if len(text) <= max_chars:
        return text
    return text[-max_chars:]


def _run_ads_api_subprocess(
    image_path: Path,
    work_dir: Path,
    model_name: str,
    model_cache_dir: Path | None = None,
    ads_python: str | Path | None = None,
    ads_package_dir: str | Path | None = None,
    ads_model_path: str | Path | None = None,
    gpu_id: int = -1,
    allow_large_images: bool = False,
) -> None:
    """Run AxonDeepSeg through its Python API in a subprocess."""
    python_exe = _resolve_ads_python(ads_python)
    env = _build_ads_env(ads_package_dir)
    model_cache_dir = model_cache_dir or (work_dir / "ads_models")
    env.update(
        {
            "TEMPTATION_ADS_IMAGE": str(image_path),
            "TEMPTATION_ADS_MODEL_NAME": model_name,
            "TEMPTATION_ADS_MODEL_DIR": str(model_cache_dir),
            "TEMPTATION_ADS_GPU_ID": str(gpu_id),
            "TEMPTATION_ADS_ALLOW_LARGE": "1" if allow_large_images else "0",
        }
    )
    if ads_model_path:
        env["TEMPTATION_ADS_MODEL_PATH"] = str(Path(ads_model_path).expanduser())

    script = r"""
import os
import inspect
from pathlib import Path

image_path = Path(os.environ["TEMPTATION_ADS_IMAGE"])
model_name = os.environ["TEMPTATION_ADS_MODEL_NAME"]
model_dir = Path(os.environ["TEMPTATION_ADS_MODEL_DIR"])
gpu_id = int(os.environ["TEMPTATION_ADS_GPU_ID"])
allow_large = os.environ.get("TEMPTATION_ADS_ALLOW_LARGE") == "1"
model_path_env = os.environ.get("TEMPTATION_ADS_MODEL_PATH")

from AxonDeepSeg.download_model import download_model
from AxonDeepSeg.segment import segment_images

if model_path_env:
    model_path = Path(model_path_env)
else:
    try:
        model_path = download_model(model_name, destination=None, overwrite=False)
    except BaseException:
        model_dir.mkdir(parents=True, exist_ok=True)
        model_path = download_model(model_name, destination=model_dir, overwrite=False)

kwargs = {
    "path_images": [image_path],
    "path_model": model_path,
    "gpu_id": gpu_id,
    "verbosity_level": 0,
}
if "allow_large_images" in inspect.signature(segment_images).parameters:
    kwargs["allow_large_images"] = allow_large
segment_images(**kwargs)
"""

    proc = subprocess.run(
        [python_exe, "-c", script],
        cwd=str(work_dir),
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AxonDeepSegError(
            "AxonDeepSeg failed.\n"
            f"Python: {python_exe}\n"
            f"Model: {model_name}\n"
            f"STDOUT:\n{_tail(proc.stdout)}\n\n"
            f"STDERR:\n{_tail(proc.stderr)}\n\n"
            "Install AxonDeepSeg in the current environment, set "
            "AXONDEEPSEG_PYTHON, or pass --ads-python pointing to an "
            "environment where AxonDeepSeg can run inference."
        )


def _first_existing(work_dir: Path, suffix: str) -> Path | None:
    matches = sorted(work_dir.glob(f"*{suffix}"))
    return matches[0] if matches else None


def _first_matching(work_dir: Path, patterns: tuple[str, ...]) -> Path | None:
    for pattern in patterns:
        matches = sorted(work_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def _read_binary_mask(path: Path | None, shape: tuple[int, int] | None = None) -> np.ndarray:
    if path is None:
        if shape is None:
            raise AxonDeepSegError("Cannot create an empty class mask before shape is known.")
        return np.zeros(shape, dtype=bool)

    arr = np.array(Image.open(path))
    if arr.ndim > 2:
        arr = arr[..., 0]
    return arr > 0


def convert_ads_masks_to_temptation_mask(
    ads_output_dir: str | Path,
    output_mask_path: str | Path,
    myelin_val: int = 64,
    axoplasm_val: int = 192,
) -> tuple[Path, dict[str, Path]]:
    """
    Convert AxonDeepSeg class masks to TEMptation's label image convention.

    AxonDeepSeg outputs binary class masks:
      *_seg-axon.png, *_seg-myelin.png, and for some models *_seg-uaxon.png.

    TEMptation expects one TIFF mask:
      64 myelin, 128 mitochondria, 192 axoplasm, 0 background.
    Mitochondria are not available from AxonDeepSeg and are left as background.
    """
    ads_output_dir = Path(ads_output_dir)
    output_mask_path = Path(output_mask_path)

    class_paths = {
        "axon": _first_existing(ads_output_dir, "_seg-axon.png"),
        "myelin": _first_existing(ads_output_dir, "_seg-myelin.png"),
        "uaxon": _first_existing(ads_output_dir, "_seg-uaxon.png")
        or _first_matching(
            ads_output_dir,
            ("*_seg-*uaxon*.png", "*_seg-*unmyelinated*axon*.png"),
        ),
    }

    if class_paths["axon"] is None and class_paths["uaxon"] is None:
        found = ", ".join(p.name for p in sorted(ads_output_dir.glob("*_seg-*.png")))
        raise AxonDeepSegError(
            "AxonDeepSeg did not produce an axon or unmyelinated-axon mask. "
            f"Found: {found or 'no *_seg-*.png files'}"
        )

    reference_path = class_paths["axon"] or class_paths["uaxon"] or class_paths["myelin"]
    reference = _read_binary_mask(reference_path)
    axon = _read_binary_mask(class_paths["axon"], shape=reference.shape)
    uaxon = _read_binary_mask(class_paths["uaxon"], shape=reference.shape)
    myelin = _read_binary_mask(class_paths["myelin"], shape=reference.shape)

    axoplasm = axon | uaxon
    out = np.zeros(axoplasm.shape, dtype=np.uint8)
    out[myelin] = np.uint8(myelin_val)
    out[axoplasm] = np.uint8(axoplasm_val)

    output_mask_path.parent.mkdir(parents=True, exist_ok=True)
    tiff.imwrite(output_mask_path, out)

    existing_class_paths = {k: v for k, v in class_paths.items() if v is not None}
    if reference_path is not None and not existing_class_paths:
        existing_class_paths["reference"] = reference_path
    return output_mask_path, existing_class_paths


def auto_segment_tem(
    tem_path: str | Path,
    output_dir: str | Path | None = None,
    model_name: str = "generalist",
    ads_python: str | Path | None = None,
    ads_package_dir: str | Path | None = None,
    ads_model_path: str | Path | None = None,
    gpu_id: int = -1,
    allow_large_images: bool = False,
    keep_intermediate: bool = False,
    output_mask_path: str | Path | None = None,
    copy_tem: bool = True,
    myelin_val: int = 64,
    axoplasm_val: int = 192,
) -> AutoSegmentationResult:
    """Run AxonDeepSeg and save a TEMptation-compatible TIFF mask."""
    tem_path = Path(tem_path).expanduser()
    if not tem_path.is_file():
        raise FileNotFoundError(f"TEM image not found: {tem_path}")
    if model_name not in ADS_MODEL_CHOICES:
        raise ValueError(
            f"model_name must be one of {', '.join(ADS_MODEL_CHOICES)}, got {model_name!r}"
        )

    output_dir = Path(output_dir).expanduser() if output_dir else tem_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    if output_mask_path is None:
        output_mask_path = output_dir / f"{tem_path.stem}_ads_mask.tif"
    else:
        output_mask_path = Path(output_mask_path).expanduser()

    output_image_path = output_dir / tem_path.name
    if copy_tem and output_image_path.resolve() != tem_path.resolve():
        shutil.copy2(tem_path, output_image_path)
    else:
        output_image_path = tem_path

    if keep_intermediate:
        work_dir = output_dir / f"{tem_path.stem}_ads_intermediate"
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True)
        cleanup = None
    else:
        cleanup = tempfile.TemporaryDirectory(prefix=f"{tem_path.stem}_ads_")
        work_dir = Path(cleanup.name)

    try:
        work_image = work_dir / tem_path.name
        shutil.copy2(tem_path, work_image)

        _run_ads_api_subprocess(
            image_path=work_image,
            work_dir=work_dir,
            model_name=model_name,
            model_cache_dir=output_dir / "ads_models",
            ads_python=ads_python,
            ads_package_dir=ads_package_dir,
            ads_model_path=ads_model_path,
            gpu_id=gpu_id,
            allow_large_images=allow_large_images,
        )

        mask_path, class_paths = convert_ads_masks_to_temptation_mask(
            ads_output_dir=work_dir,
            output_mask_path=output_mask_path,
            myelin_val=myelin_val,
            axoplasm_val=axoplasm_val,
        )

        result_work_dir = work_dir if keep_intermediate else None
        return AutoSegmentationResult(
            image_path=output_image_path,
            mask_path=mask_path,
            model_name=model_name,
            class_mask_paths=class_paths,
            output_dir=output_dir,
            ads_work_dir=result_work_dir,
        )
    finally:
        if cleanup is not None:
            cleanup.cleanup()
