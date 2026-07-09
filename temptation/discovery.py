"""Find TEM/mask pairs and infer groups.

Phase 1: `find_pairs` is a verbatim move of the original
`_find_pairs_in_folder` (non-recursive, *.tif/*.tiff only, warnings.warn on
skip). Recursive discovery, the extended glob (covers .gif etc.), and
group inference (F7 fixes from IMPLEMENTATION_BLUEPRINT.md Sec 0) land in
Phase 2 -- do not add them here without also updating the Phase 1 golden
regression tests, which assert this exact (limited) behavior.
"""

import re
import warnings
from pathlib import Path


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
