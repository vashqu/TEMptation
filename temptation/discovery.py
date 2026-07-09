"""Find TEM/mask pairs and infer groups.

`find_pairs` is a verbatim move of the original `_find_pairs_in_folder`
(non-recursive, *.tif/*.tiff only, warnings.warn on skip) and is left
untouched -- it is still what plain `--folder X` uses, and the Phase 1
golden regression tests assert this exact (limited) behavior.

Phase 2 adds a *separate* code path for the F7 fixes (recursive scan,
extended glob covering .gif/.png, skips reported as data instead of
warnings): `find_groups` and its helper `_scan_folder_recursive`. This
powers the new --groups/--input-root/--recursive CLI flags without
changing what `--folder X` alone does.
"""

import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


IMAGE_GLOB_PATTERNS = ("*.tif", "*.tiff", "*.png", "*.gif")


@dataclass
class SkippedFile:
    path: Path
    reason: str


def find_pairs(folder: Path):
    """
    Scan folder for paired TEM + mask files.

    Supported naming conventions
    ----------------------------
    Pattern A -- explicit prefix with underscore separator:
        tem_<id>.tif   OR   axon_<id>.tif   ->  TEM image
        mask_<id>.tif                        ->  mask image

    Pattern B -- numeric prefix with tissue keyword anywhere in name:
        <number>. Axon *.tif  OR  <number>. TEM *.tif  ->  TEM image
        <number>. Mask *.tif                            ->  mask image
        IDs are matched by the leading integer.

    Returns list of dicts: {id, tem_path, mask_path}
    """
    tif_files = sorted(folder.glob("*.tif")) + sorted(folder.glob("*.tiff"))

    tem_map = {}
    mask_map = {}

    for f in tif_files:
        stem = f.stem
        stem_lower = stem.lower()

        if stem_lower.startswith("tem_") or stem_lower.startswith("axon_"):
            img_id = stem.split("_", 1)[1]
            tem_map.setdefault(img_id, f)
        elif stem_lower.startswith("mask_"):
            img_id = stem.split("_", 1)[1]
            mask_map.setdefault(img_id, f)
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
            "id": img_id,
            "tem_path": tem_map[img_id],
            "mask_path": mask_map[img_id],
        })

    return pairs


def _scan_folder_recursive(folder: Path):
    """Recursive scan of `folder` (and subdirectories) for TEM+mask pairs,
    using the extended glob (tif/tiff/png/gif -- fixes F7's silently
    dropped .gif masks). Returns (pairs: list[dict], skipped: list[SkippedFile])
    -- skips are returned as data rather than emitted via warnings.warn,
    since this path is meant for unattended batch/group runs where an
    interactive warning is easy to miss."""
    image_files = []
    seen = set()
    for pattern in IMAGE_GLOB_PATTERNS:
        for f in sorted(folder.rglob(pattern)):
            if f not in seen:
                seen.add(f)
                image_files.append(f)

    tem_map = {}
    mask_map = {}
    skipped = []

    for f in image_files:
        stem = f.stem
        stem_lower = stem.lower()

        if stem_lower.startswith("tem_") or stem_lower.startswith("axon_"):
            img_id = stem.split("_", 1)[1]
            tem_map.setdefault(img_id, f)
        elif stem_lower.startswith("mask_"):
            img_id = stem.split("_", 1)[1]
            mask_map.setdefault(img_id, f)
        else:
            m = re.match(r"^(\d+)", stem)
            if not m:
                skipped.append(SkippedFile(path=f, reason="no leading id and no tem_/axon_/mask_ prefix"))
                continue
            img_id = m.group(1)
            if "mask" in stem_lower:
                mask_map.setdefault(img_id, f)
            elif "axon" in stem_lower or "tem" in stem_lower:
                tem_map.setdefault(img_id, f)
            else:
                skipped.append(SkippedFile(path=f, reason="leading id but no axon/tem/mask keyword in name"))

    def _sort_key(x):
        return int(x) if x.isdigit() else x

    pairs = []
    all_ids = set(tem_map) | set(mask_map)
    for img_id in sorted(all_ids, key=_sort_key):
        if img_id not in tem_map:
            skipped.append(SkippedFile(path=mask_map[img_id], reason=f"no TEM file for id={img_id!r}"))
            continue
        if img_id not in mask_map:
            skipped.append(SkippedFile(path=tem_map[img_id], reason=f"no mask file for id={img_id!r}"))
            continue
        pairs.append({
            "id": img_id,
            "tem_path": tem_map[img_id],
            "mask_path": mask_map[img_id],
        })

    return pairs, skipped


def find_groups(input_root: Path, group_map: dict):
    """Scan input_root/<folder_name> recursively for each entry in
    group_map ({folder_name: group_label}), tagging every discovered pair
    with its group. Returns (pairs: list[dict], skipped: list[SkippedFile]).

    Example: find_groups(Path("."), {"normal_data": "normal",
    "pathological_data": "pathological"})
    """
    input_root = Path(input_root)
    all_pairs = []
    all_skipped = []

    for folder_name, group_label in group_map.items():
        folder = input_root / folder_name
        if not folder.is_dir():
            all_skipped.append(SkippedFile(path=folder, reason="group folder not found"))
            continue
        pairs, skipped = _scan_folder_recursive(folder)
        for p in pairs:
            p["group"] = group_label
            p["source_folder"] = folder
        all_pairs.extend(pairs)
        all_skipped.extend(skipped)

    return all_pairs, all_skipped


def infer_group(path: Path, group_map: dict) -> Optional[str]:
    """Walk path's parts looking for a folder name present in group_map.
    Returns the mapped group label, or None if no part matches."""
    for part in Path(path).parts:
        if part in group_map:
            return group_map[part]
    return None


def parse_group_map(specs) -> dict:
    """Parse ["normal_data:normal", "pathological_data:pathological"] into
    {"normal_data": "normal", "pathological_data": "pathological"}."""
    group_map = {}
    for spec in specs:
        if ":" not in spec:
            raise ValueError(f"--groups entries must be FOLDER:LABEL, got {spec!r}")
        folder_name, label = spec.split(":", 1)
        group_map[folder_name] = label
    return group_map
