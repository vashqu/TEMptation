"""Export panel (blueprint Sec 8.8): a literal checklist of what will
be written, with resolved absolute paths, row-count estimates, and an
overwrite warning per file, plus a single Export button. This panel
assembles a completed run's already-computed DataFrames into files by
calling the same temptation.export writers the CLI uses -- it computes
nothing itself.

Plots/ (per-image overlay PNGs) is deliberately not offered here:
analyze_dataset doesn't retain each image's labels_ws (a memory
tradeoff -- see pipeline.py), so overlay export needs a per-image
re-analysis, which belongs with the visual review panel (Phase 7d),
not this one."""

import dataclasses
import tkinter as tk
from tkinter import messagebox, ttk

from temptation import export
from temptation import qc as qc_mod
from temptation import summaries

from ..widgets import CARD_BG, LBL_FG, SUCCESS_FG, WARN_FG, make_card


def build(parent, state):
    parent.columnconfigure(0, weight=1)

    outer, inner = make_card(parent, "Export")
    outer.grid(row=0, column=0, sticky="ew")
    inner.columnconfigure(0, weight=1)

    status_var = tk.StringVar(value="")
    tk.Label(inner, textvariable=status_var, bg=CARD_BG, fg=LBL_FG,
             justify="left", anchor="w", wraplength=680).pack(anchor="w", pady=(0, 10))

    rows_frame = tk.Frame(inner, bg=CARD_BG)
    rows_frame.pack(fill="x")

    row_vars: dict[str, tk.BooleanVar] = {}

    export_btn = tk.Button(inner, text="Export files", state="disabled")
    export_btn.pack(anchor="w", pady=(12, 0))

    def _refresh():
        for w in rows_frame.winfo_children():
            w.destroy()
        row_vars.clear()

        result = state.result
        if result is None or state.output_dir is None:
            status_var.set("Run analysis first (⑤ Run), with an output directory set (① Load data).")
            export_btn.configure(state="disabled", command=lambda: None)
            return

        out_dir = state.output_dir
        specs = _row_specs(state, result, out_dir)

        for i, spec in enumerate(specs):
            var = tk.BooleanVar(value=spec["default_on"])
            row_vars[spec["key"]] = var
            cb = ttk.Checkbutton(rows_frame, variable=var, style="Card.TCheckbutton")
            if spec["locked"]:
                cb.configure(state="disabled")
            cb.grid(row=i, column=0, sticky="w", pady=2)

            tk.Label(rows_frame, text=spec["label"], bg=CARD_BG, fg=LBL_FG,
                     font=("TkDefaultFont", 9, "bold"), width=22, anchor="w").grid(
                row=i, column=1, sticky="w")

            count_text = f"~{spec['n_rows']} row(s)" if spec["n_rows"] is not None else ""
            tk.Label(rows_frame, text=count_text, bg=CARD_BG, fg="#999999", width=12,
                     font=("TkDefaultFont", 9)).grid(row=i, column=2, sticky="w", padx=(10, 10))

            note_color = WARN_FG if spec["warn"] else ("#999999" if spec["note"] else SUCCESS_FG)
            tk.Label(rows_frame, text=spec["note"], bg=CARD_BG, fg=note_color,
                     font=("TkDefaultFont", 8), width=34, anchor="w").grid(row=i, column=3, sticky="w")

            tk.Label(rows_frame, text=str(spec["path"]), bg=CARD_BG, fg="#BBBBBB",
                     font=("TkDefaultFont", 8)).grid(row=i, column=4, sticky="w", padx=(10, 0))

        status_var.set(f"Output directory: {out_dir}")
        export_btn.configure(
            state="normal",
            command=lambda: _do_export(state, result, out_dir, specs, row_vars, status_var),
        )

    state.on_change(_refresh)
    _refresh()


def _row_specs(state, result, out_dir):
    """One dict per candidate output file. `locked` rows can't be
    unchecked (axons/image are always written; the manifest is always
    written and isn't even offered as a checkbox row for that reason)."""
    n_groups = len({p.get("group") for p in state.all_pairs() if p.get("group")})

    specs = [
        dict(key="axons", label="axons.csv", path=out_dir / "axons.csv",
             n_rows=len(result.df_axons), default_on=True, locked=True,
             note="always written", warn=False),
        dict(key="image", label="image_summary.csv", path=out_dir / "image_summary.csv",
             n_rows=len(result.df_image), default_on=True, locked=True,
             note="always written", warn=False),
    ]

    if state.collect_mito and not result.df_mito.empty:
        specs.append(dict(
            key="mito", label="mitochondria_metrics.csv", path=out_dir / "mitochondria_metrics.csv",
            n_rows=len(result.df_mito), default_on=True, locked=False, note="", warn=False,
        ))
    else:
        specs.append(dict(
            key="mito", label="mitochondria_metrics.csv", path=out_dir / "mitochondria_metrics.csv",
            n_rows=None, default_on=False, locked=True,
            note="enable in ③ Metrics ('Collect mitochondria_metrics.csv')", warn=False,
        ))

    qc_df = qc_mod.qc_report(result.df_axons) if not result.df_axons.empty else result.df_qc_report
    specs.append(dict(
        key="qc", label="qc_report.csv", path=out_dir / "qc_report.csv",
        n_rows=len(qc_df), default_on=True, locked=False, note="", warn=False,
    ))

    if n_groups >= 1:
        group_df = summaries.summarize_groups(result.df_image)
        specs.append(dict(
            key="group", label="group_metrics.csv", path=out_dir / "group_metrics.csv",
            n_rows=len(group_df), default_on=n_groups > 1, locked=False,
            note="" if n_groups > 1 else "only one group loaded", warn=False,
        ))
    else:
        specs.append(dict(
            key="group", label="group_metrics.csv", path=out_dir / "group_metrics.csv",
            n_rows=None, default_on=False, locked=True, note="no group labels present", warn=False,
        ))

    for spec in specs:
        if spec["path"].exists():
            spec["warn"] = True
            spec["note"] = "⚠ will overwrite existing file" if not spec["note"] else spec["note"] + " -- ⚠ exists"

    return specs


def _do_export(state, result, out_dir, specs, row_vars, status_var):
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    row_counts = {"axons": len(result.df_axons), "images": len(result.df_image),
                  "mitochondria": len(result.df_mito), "qc_report": len(result.df_qc_report),
                  "groups": 0}

    try:
        if row_vars["axons"].get():
            outputs["axons_csv"] = export.write_axon_csv(result.df_axons, out_dir)
        if row_vars["image"].get():
            outputs["image_summary_csv"] = export.write_image_csv(result.df_image, out_dir)
        if row_vars["mito"].get() and not result.df_mito.empty:
            outputs["mitochondria_metrics_csv"] = export.write_mito_csv(result.df_mito, out_dir)
        if row_vars["qc"].get():
            qc_df = qc_mod.qc_report(result.df_axons) if not result.df_axons.empty else result.df_qc_report
            outputs["qc_report_csv"] = export.write_qc_report_csv(qc_df, out_dir)
        if row_vars["group"].get():
            group_df = summaries.summarize_groups(result.df_image)
            row_counts["groups"] = len(group_df)
            outputs["group_metrics_csv"] = export.write_group_csv(group_df, out_dir)

        manifest_config = {
            "segmentation": dataclasses.asdict(state.seg_cfg),
            "qc_thresholds": dataclasses.asdict(state.qc_thresholds),
            "exclude_qc_failed": state.exclude_qc_failed,
            "collect_mito": state.collect_mito,
            "pixel_size_um": state.pixel_size_um,
            "pixel_space_only": state.pixel_space_only,
            "resolved_output_dir": str(out_dir),
        }
        manifest_inputs = [
            {"id": p["id"], "tem_path": str(p["tem_path"]), "mask_path": str(p["mask_path"]),
             "group": p.get("group")}
            for p in state.all_pairs()
        ]
        manifest_path = export.write_run_manifest(
            out_dir,
            config=manifest_config,
            inputs=manifest_inputs,
            skipped_files=[
                {"path": str(s.path), "reason": s.reason}
                for entry in state.groups.values() for s in entry.skipped
            ],
            outputs=outputs,
            row_counts=row_counts,
            exclusion_counts=qc_mod.exclusion_reason_counts(result.df_axons),
            processing_errors=[{"id": pid, "message": msg} for pid, msg, _tb in result.errors],
        )
        outputs["run_manifest_json"] = manifest_path

        state.step_status["export"] = "complete"
        # notify() first (rebuilds the checklist -- overwrite warnings
        # clear once a file has just been (re)written), then set the
        # success message so the rebuild doesn't stomp it.
        state.notify()
        status_var.set("Exported:\n" + "\n".join(f"  {k} → {v}" for k, v in outputs.items()))
    except Exception as exc:
        messagebox.showerror("Export failed", str(exc))
