"""CLI: argparse, batch workflow, main(). Phase 1 behavior (--folder /
--tem+--mask) is untouched; Phase 2 adds --input-root/--groups/--group/
--recursive for multi-group batch runs with group labels and metadata
columns (IMPLEMENTATION_BLUEPRINT.md Sec 5.4/5.6, fixes F3/F7)."""

import argparse
import sys
import traceback
from pathlib import Path

import pandas as pd

from . import dataio
from . import discovery
from . import export
from . import pipeline
from . import plotting
from . import schema
from .config import SegmentationConfig


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Segment and measure nerve fibers from TEM + mask images.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    input_grp = p.add_argument_group(
        "Input (use --folder OR --tem + --mask OR --groups)"
    )
    input_grp.add_argument("--folder", type=Path, default=None)
    input_grp.add_argument("--tem", type=Path, default=None)
    input_grp.add_argument("--mask", type=Path, default=None)
    input_grp.add_argument(
        "--input-root", type=Path, default=None,
        help="Root directory that --groups folder names are resolved relative to "
             "(defaults to current directory).",
    )
    input_grp.add_argument(
        "--groups", nargs="+", default=None, metavar="FOLDER:LABEL",
        help="Batch-process multiple labeled groups in one run, e.g. "
             "--groups normal_data:normal pathological_data:pathological. "
             "Each folder is scanned recursively.",
    )
    input_grp.add_argument(
        "--group", type=str, default=None,
        help="Explicit group label to stamp on every row (used with --folder "
             "or --tem/--mask; ignored with --groups, which assigns groups "
             "per folder).",
    )
    input_grp.add_argument(
        "--recursive", action="store_true",
        help="With --folder, scan subdirectories too and use the extended "
             "glob (tif/tiff/png/gif) instead of the default non-recursive "
             "tif/tiff-only scan.",
    )

    scale_grp = p.add_argument_group("Scale (provide exactly one option)")
    scale_ex = scale_grp.add_mutually_exclusive_group()
    scale_ex.add_argument("--pixel-um", type=float, default=None)
    scale_ex.add_argument("--bar", nargs=2, metavar=("PX", "UM"), type=float, default=None)

    p.add_argument("--mode", choices=["normal", "pathological", "auto"], default="auto")
    p.add_argument("--myelin-threshold", type=int, default=200)
    p.add_argument("--myelin-val", type=int, default=64)
    p.add_argument("--axoplasm-val", type=int, default=192)
    p.add_argument("--mito-val", type=int, default=128)
    p.add_argument("--smoothing-radius", type=int, default=1)
    p.add_argument("--min-axon-area", type=int, default=200)
    p.add_argument("--min-myelin-area", type=int, default=300)
    p.add_argument("--watershed-mode", choices=["simple", "weighted"], default="weighted")
    p.add_argument("--watershed-weight", choices=["radius", "area"], default="radius")
    p.add_argument("--watershed-compactness", type=float, default=0.001)
    p.add_argument(
        "--watershed-beta", type=float, default=1.0,
        help="Bias strength for weighted watershed; higher = larger axons claim more territory.",
    )
    p.add_argument(
        "--assign-detached-myelin", choices=["none", "nearest"], default="none",
        help="Assign detached myelin pixels to the nearest axon for myelin metric computation.",
    )
    p.add_argument(
        "--mito-hole-handling", choices=["fill", "legacy"], default="fill",
        help="'fill' (default) correctly attributes a full mitochondrion to its "
             "axon by connectivity; 'legacy' reproduces the pre-fix behavior "
             "that only captured a ~1px rim of each mitochondrion (see "
             "IMPLEMENTATION_BLUEPRINT.md Sec 0, defect F1) -- use 'legacy' "
             "only to reproduce numbers from before this fix.",
    )
    p.add_argument(
        "--mito-assignment", choices=["centroid", "overlap", "legacy"], default="centroid",
        help="How each mitochondrion is attributed to an axon. 'centroid' "
             "(default) labels mitochondria once globally and assigns each to "
             "the axon at its centroid; 'overlap' uses maximum pixel overlap "
             "instead; 'legacy' reproduces the pre-fix per-fiber-crop "
             "intersection, which fragments and double-counts any "
             "mitochondrion straddling a watershed boundary (defect F8) -- "
             "use 'legacy' only to reproduce numbers from before this fix.",
    )
    p.add_argument(
        "--demyelination-reference", type=float, default=None, metavar="FRACTION",
        help="Reference myelin_area_fraction_of_fov used to compute "
             "image_demyelination_index = 1 - (this image's fraction / reference), "
             "clipped to [0,1]. No default: this is a biological calibration choice "
             "(IMPLEMENTATION_BLUEPRINT.md Sec 11 A4) -- image_demyelination_index "
             "is NaN unless you supply this explicitly, e.g. from the median "
             "myelin_area_fraction_of_fov of a known-normal reference set.",
    )
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--plot", action="store_true")
    p.add_argument(
        "--write-mito-csv", action="store_true",
        help="Also write mitochondria_metrics.csv (one row per real mitochondrion). "
             "Requires --mito-assignment centroid or overlap; ignored (with a warning) "
             "under --mito-assignment legacy, since a per-fiber-crop mitochondrion "
             "fragment is not a real, whole mitochondrion.",
    )

    return p


def resolve_pixel_length(args) -> float:
    if args.pixel_um is not None:
        return args.pixel_um
    elif args.bar is not None:
        bar_px, bar_um = args.bar
        return bar_um / bar_px
    else:
        raise SystemExit(
            "ERROR: You must provide either --pixel-um or --bar <PX> <UM> to set the scale."
        )


def process_pair(
    image_id: str,
    tem_path: Path,
    mask_path: Path,
    pixel_length_um: float,
    args,
    output_dir: Path,
    group: str = None,
) -> tuple:
    tem = dataio.read_image(tem_path)
    mask = dataio.read_mask(mask_path)

    seg_cfg = SegmentationConfig(
        myelin_val=args.myelin_val,
        axoplasm_val=args.axoplasm_val,
        mito_val=args.mito_val,
        mode=args.mode,
        myelin_threshold_px=args.myelin_threshold,
        smoothing_radius_px=args.smoothing_radius,
        min_axon_area_px=args.min_axon_area,
        min_myelin_area_px=args.min_myelin_area,
        watershed_mode=args.watershed_mode,
        watershed_weight=args.watershed_weight,
        watershed_compactness=args.watershed_compactness,
        watershed_beta=args.watershed_beta,
        assign_detached_myelin=args.assign_detached_myelin,
        mito_hole_handling=args.mito_hole_handling,
        mito_assignment=args.mito_assignment,
    )
    df_axons, df_image, labels_ws, resolved_mode, df_mito = pipeline.analyze_image_legacy(
        tem, mask, pixel_length_um, seg_cfg,
        reference_myelin_fraction=args.demyelination_reference,
    )

    print(f"  [{resolved_mode.upper()}] id={image_id!r}  "
          f"→ {len(df_axons)} axons detected")

    for df in (df_axons, df_image):
        if not df.empty:
            if "image_id" not in df.columns:
                df.insert(0, "image_id", image_id)
            else:
                df["image_id"] = image_id
            if "mode" not in df.columns:
                df.insert(1, "mode", resolved_mode)
            else:
                df["mode"] = resolved_mode
            # Identity/metadata columns (Phase 2). Appended, not inserted,
            # so the legacy column positions above are undisturbed -- see
            # schema.py's note on why the full identity-block-prepend
            # redesign is deferred to Phase 6.
            df["group"] = group
            df["image_path"] = str(tem_path)
            df["mask_path"] = str(mask_path)
            df["pixel_size_um"] = pixel_length_um
            df["schema_version"] = schema.resolve_schema_version(args.mito_hole_handling)

    if getattr(args, "write_mito_csv", False):
        if not df_mito.empty:
            df_mito.insert(0, "image_id", image_id)
            df_mito["group"] = group
        elif args.mito_assignment == "legacy":
            print(f"  [WARNING] --write-mito-csv has no effect under "
                  f"--mito-assignment legacy for id={image_id!r} "
                  f"(no per-mitochondrion table in that mode).")

    if args.plot and not df_axons.empty:
        plot_path = output_dir / f"overlay_{image_id}.png"
        try:
            fig = plotting.overlay_figure(
                tem=tem,
                mask=mask,
                labels_ws=labels_ws,
                df_axons=df_axons,
                myelin_val=args.myelin_val,
                axoplasm_val=args.axoplasm_val,
                mito_val=args.mito_val,
                title=f"Image {image_id}  [{resolved_mode}]",
            )
            plotting.save_overlay(fig, plot_path)
            plotting.show_overlay_nonblocking(fig)
        except Exception as plot_exc:
            print(f"  [WARNING] Plot failed for id={image_id!r}: {plot_exc}")
            traceback.print_exc()

    return df_axons, df_image, df_mito


def main():
    parser = build_parser()
    args = parser.parse_args()

    skipped_files = []

    if args.groups is not None:
        if args.folder is not None or args.tem is not None or args.mask is not None:
            parser.error("Use --groups OR --folder OR (--tem + --mask), not together.")
        if args.group is not None:
            parser.error(
                "--group is ignored with --groups (each folder in --groups "
                "already carries its own label); remove one of the two."
            )
        group_map = discovery.parse_group_map(args.groups)
        input_root = args.input_root or Path(".")
        if not input_root.is_dir():
            parser.error(f"--input-root {input_root!r} is not a directory.")
        pairs, skipped_files = discovery.find_groups(input_root, group_map)
        if not pairs:
            sys.exit(f"ERROR: No valid TEM+mask pairs found under {input_root!r} for groups {group_map}")
        output_dir = args.output_dir or input_root
    elif args.folder is not None:
        if args.tem is not None or args.mask is not None:
            parser.error("Use either --folder OR (--tem + --mask), not both.")
        if not args.folder.is_dir():
            parser.error(f"--folder {args.folder!r} is not a directory.")
        if args.recursive:
            pairs, skipped_files = discovery._scan_folder_recursive(args.folder)
        else:
            pairs = discovery.find_pairs(args.folder)
        if not pairs:
            sys.exit(f"ERROR: No valid TEM+mask pairs found in {args.folder}")
        for pair in pairs:
            pair.setdefault("group", args.group)
        output_dir = args.output_dir or args.folder
    else:
        if args.tem is None or args.mask is None:
            parser.error("Provide --folder OR both --tem and --mask OR --groups.")
        if not args.tem.is_file():
            parser.error(f"--tem file not found: {args.tem}")
        if not args.mask.is_file():
            parser.error(f"--mask file not found: {args.mask}")
        pairs = [{"id": args.tem.stem, "tem_path": args.tem, "mask_path": args.mask, "group": args.group}]
        output_dir = args.output_dir or args.tem.parent

    output_dir.mkdir(parents=True, exist_ok=True)
    pixel_length_um = resolve_pixel_length(args)
    print(f"Pixel length: {pixel_length_um:.6f} µm/px")

    if skipped_files:
        print(f"\n[WARNING] {len(skipped_files)} file(s) skipped during discovery:")
        for s in skipped_files:
            print(f"  SKIP: {s.path}  ({s.reason})")

    all_axons = []
    all_images = []
    all_mito = []

    for pair in pairs:
        try:
            df_axons, df_image, df_mito = process_pair(
                image_id=pair["id"],
                tem_path=pair["tem_path"],
                mask_path=pair["mask_path"],
                pixel_length_um=pixel_length_um,
                args=args,
                output_dir=output_dir,
                group=pair.get("group"),
            )
            if not df_axons.empty:
                all_axons.append(df_axons)
            if not df_image.empty:
                all_images.append(df_image)
            if not df_mito.empty:
                all_mito.append(df_mito)
        except Exception as exc:
            print(f"\n[ERROR] id={pair['id']!r} failed with: {exc}")
            traceback.print_exc()
            print()

    if all_axons:
        df_all_axons = pd.concat(all_axons, ignore_index=True)
        axon_csv = export.write_axon_csv(df_all_axons, output_dir)
        print(f"\nAxon-level results  → {axon_csv}  ({len(df_all_axons)} rows)")
    else:
        print("\nNo axons detected across all images.")
        df_all_axons = pd.DataFrame()

    if all_images:
        df_all_images = pd.concat(all_images, ignore_index=True)
        img_csv = export.write_image_csv(df_all_images, output_dir)
        print(f"Image-level summary → {img_csv}  ({len(df_all_images)} rows)")

    if args.write_mito_csv:
        if all_mito:
            df_all_mito = pd.concat(all_mito, ignore_index=True)
            mito_csv = export.write_mito_csv(df_all_mito, output_dir)
            print(f"Mitochondria results → {mito_csv}  ({len(df_all_mito)} rows)")
        else:
            print("No mitochondria to write (empty or --mito-assignment legacy).")

    return df_all_axons
