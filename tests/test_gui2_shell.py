"""Smoke tests for the Phase 7 GUI shell (gui/app.py): construct the
real Tk app, drive its dataset/calibration state through the same
methods the widgets call, and confirm Run wiring reaches
temptation.pipeline.analyze_dataset with results structurally
equivalent to the CLI path. Skips gracefully with no display."""

import os
import time
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

DATA_ROOT = Path(__file__).parent.parent.parent
NORMAL_DIR = DATA_ROOT / "normal_data" / "162-165"


@pytest.fixture
def app():
    import tkinter as tk

    from gui.app import TEMptationApp

    try:
        instance = TEMptationApp()
    except tk.TclError as exc:
        pytest.skip(f"no display available for Tk: {exc}")
    yield instance
    instance.destroy()


def test_gui2_constructs_with_all_steps(app):
    assert app.title() == "TEMptation — Nerve Morphometry"
    assert set(app._panel_frames.keys()) == {
        "dataset", "calibration", "metrics", "qc", "run", "review", "export",
    }


def test_gui2_step_navigation_updates_rail(app):
    app._select_step("calibration")
    assert app._current_step == "calibration"
    app._select_step("run")
    assert app._current_step == "run"


def test_gui2_dataset_scan_populates_state(app, tmp_path):
    from temptation import discovery

    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, skipped = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.folder = NORMAL_DIR
    entry.pairs = pairs
    entry.skipped = skipped
    app.state_.groups["normal"] = entry
    app.state_.output_dir = tmp_path
    app.state_.pixel_size_um = 0.00524
    app.state_.mark_inputs_changed()

    assert app.state_.total_pairs() == 4
    assert app.state_.step_status["run"] == "not_started"


def test_gui2_run_calls_analyze_dataset_and_matches_cli(app, tmp_path):
    """The one that matters (blueprint Sec 10.6): drive a real run
    through the GUI's own _on_run/_run_worker path and confirm the
    resulting axon count matches what the CLI produces for the same
    folder, proving Run doesn't reimplement or diverge from
    pipeline.analyze_dataset."""
    from temptation import discovery

    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, skipped = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.folder = NORMAL_DIR
    entry.pairs = pairs
    app.state_.groups["normal"] = entry
    app.state_.output_dir = tmp_path
    app.state_.pixel_size_um = 0.00524

    app._on_run()

    # Drive Tk's real mainloop rather than busy-polling app.update() --
    # a tight update() loop racing the worker thread's cross-thread
    # after(0, ...) calls segfaults Tcl under macOS Aqua. Scheduling a
    # periodic after() check and quitting the loop once done is the
    # standard pattern and lets Tk process those callbacks normally.
    deadline = time.time() + 30

    def _poll():
        if app.state_.run_in_progress and time.time() < deadline:
            app.after(50, _poll)
        else:
            app.quit()

    app.after(50, _poll)
    app.mainloop()

    assert not app.state_.run_in_progress
    result = app.state_.result
    assert result is not None
    assert len(result.df_axons) == 47  # matches tests/golden_v2/normal row count
    assert len(result.errors) == 0
    assert app.state_.step_status["run"] == "complete"
    # last_raw_axons must survive a subsequent mark_inputs_changed() (e.g.
    # editing a QC threshold) -- it's what the QC panel's live preview
    # recomputes against without a re-run (blueprint Sec 8.5).
    assert app.state_.last_raw_axons is not None
    app.state_.mark_inputs_changed()
    assert app.state_.last_raw_axons is not None
    assert app.state_.result is None


def test_gui2_metrics_scan_detects_myelin_presence(app):
    """normal_data/162-165's masks all contain myelin -- the scan should
    report full availability, not just "not yet scanned"."""
    from temptation import discovery

    from gui.panels.metrics import _scan
    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.pairs = pairs
    app.state_.groups["normal"] = entry

    import tkinter as tk
    availability_var = tk.StringVar()
    status_var = tk.StringVar()
    status_lbl = tk.Label(app)
    card_widgets = {"g-ratio distribution": (status_var, status_lbl, ("myelin",), False)}

    _scan(app.state_, availability_var, card_widgets)

    assert "myelin in 4" in availability_var.get()
    assert "available" in status_var.get()


def test_gui2_results_and_export_after_real_run(app, tmp_path):
    """Drive a real run through the GUI, then build the results
    dashboard (exercising the matplotlib boxplot/jitter path) and the
    export panel, and actually click Export -- confirming the checklist
    writes real files via temptation.export, not a GUI-side
    reimplementation."""
    import tkinter as tk

    from temptation import discovery

    from gui.panels import export as export_panel
    from gui.panels import results as results_panel
    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.pairs = pairs
    app.state_.groups["normal"] = entry
    app.state_.output_dir = tmp_path
    app.state_.pixel_size_um = 0.00524

    app._on_run()
    deadline = time.time() + 30

    def _poll():
        if app.state_.run_in_progress and time.time() < deadline:
            app.after(50, _poll)
        else:
            app.quit()

    app.after(50, _poll)
    app.mainloop()
    assert app.state_.result is not None

    results_panel.build(tk.Frame(app), app.state_)

    export_frame = tk.Frame(app)
    export_panel.build(export_frame, app.state_)

    def _find_button(widget):
        if isinstance(widget, tk.Button) and widget.cget("text") == "Export files":
            return widget
        for child in widget.winfo_children():
            found = _find_button(child)
            if found is not None:
                return found
        return None

    export_btn = _find_button(export_frame)
    assert export_btn is not None
    assert str(export_btn.cget("state")) == "normal"
    export_btn.invoke()

    assert (tmp_path / "axons.csv").exists()
    assert (tmp_path / "image_summary.csv").exists()
    assert (tmp_path / "run_manifest.json").exists()


def test_gui2_metrics_and_qc_panels_build_with_prior_run_data(app):
    """Building the metrics/QC panels with a populated dataset and a
    prior run's raw axons exercises the mask-scan setup and the QC
    live-preview recompute path, both of which only run once real data
    exists -- a bare construction test wouldn't reach either."""
    import tkinter as tk

    import pandas as pd
    from temptation import discovery

    from gui.panels import metrics as metrics_panel
    from gui.panels import qc as qc_panel
    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.pairs = pairs
    app.state_.groups["normal"] = entry
    app.state_.last_raw_axons = pd.read_csv(
        Path(__file__).parent / "golden_v2" / "normal" / "axons.csv"
    )

    metrics_panel.build(tk.Frame(app), app.state_)
    qc_panel.build(tk.Frame(app), app.state_)
