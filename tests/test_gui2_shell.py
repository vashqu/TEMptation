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
