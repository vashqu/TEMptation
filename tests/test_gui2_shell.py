"""Smoke tests for the TEMptation GUI shell (gui/app.py): construct the
real Tk app, drive its dataset/calibration state through the same
methods the widgets call, and confirm Run wiring reaches
temptation.pipeline.analyze_dataset with results structurally
equivalent to the CLI path. Skips gracefully with no display.

Covers the post-launch reorg (7 rail steps -> 3: Setup/Review/Export,
Setup bundling Data/Calibration/Segmentation/QC as sub-tabs with a
pinned Run footer) and the accompanying fixes: state.result surviving
input edits (result_stale instead of nulling it), debounced
calibration/QC input handlers, the Inspect panel's list-squeeze layout
bug, and the optional 'mode' column."""

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


def _run_and_wait(app, pairs, group_name, output_dir, pixel_um=0.00524):
    """Shared helper: populate one group, drive a real batch run through
    app._on_run()/_run_worker(), and block (via the after()+mainloop()/
    quit() pattern -- a tight app.update() loop is known to segfault
    Tcl under macOS Aqua when racing the worker thread's cross-thread
    after(0, ...) calls) until it completes."""
    from gui.state import GroupEntry

    entry = GroupEntry(name=group_name)
    entry.pairs = pairs
    app.state_.groups[group_name] = entry
    app.state_.output_dir = output_dir
    app.state_.pixel_size_um = pixel_um

    app._on_run()
    deadline = time.time() + 30

    def _poll():
        if app.state_.run_in_progress and time.time() < deadline:
            app.after(50, _poll)
        else:
            app.quit()

    app.after(50, _poll)
    app.mainloop()


def test_gui2_constructs_with_three_steps(app):
    assert app.title() == "TEMptation — Nerve Morphometry"
    assert set(app._panel_frames.keys()) == {"setup", "review", "export"}


def test_gui2_step_navigation_updates_rail(app):
    app._select_step("review")
    assert app._current_step == "review"
    app._select_step("export")
    assert app._current_step == "export"
    app._select_step("setup")
    assert app._current_step == "setup"


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
    assert app.state_.step_status["review"] == "not_started"


def test_gui2_run_calls_analyze_dataset_and_matches_cli(app, tmp_path):
    """The one that matters (blueprint Sec 10.6): drive a real run
    through the GUI's own _on_run/_run_worker path and confirm the
    resulting axon count matches what the CLI produces for the same
    folder, proving Run doesn't reimplement or diverge from
    pipeline.analyze_dataset."""
    from temptation import discovery

    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    _run_and_wait(app, pairs, "normal", tmp_path)

    assert not app.state_.run_in_progress
    result = app.state_.result
    assert result is not None
    assert len(result.df_axons) == 47  # matches tests/golden_v2/normal row count
    assert len(result.errors) == 0
    assert app.state_.step_status["setup"] == "complete"
    assert app.state_.result_stale is False


def test_gui2_result_survives_input_edits_but_flips_stale(app, tmp_path):
    """The behavior this whole reorg fixed: editing an input used to
    null state.result outright (blanking the dashboard/export checklist
    mid-keystroke). Now it must survive, with result_stale flipping to
    tell the user a fresh Run would reflect the new parameters."""
    from temptation import discovery

    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    _run_and_wait(app, pairs, "normal", tmp_path)
    assert app.state_.result is not None
    assert app.state_.result_stale is False

    app.state_.mark_inputs_changed()

    assert app.state_.result is not None
    assert app.state_.result_stale is True
    assert app.state_.step_status["review"] == "not_started"


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


def test_gui2_segmentation_mode_selector_updates_state(app):
    """Metrics panel's mode Combobox drives state.seg_cfg.mode directly
    -- confirms the control promoted out of the generic parameter
    accordion actually wires to the same field the backend reads."""
    import tkinter as tk
    from tkinter import ttk

    from gui.panels import metrics as metrics_panel

    frame = tk.Frame(app)
    metrics_panel.build(frame, app.state_)

    def _find(widget, cls):
        if isinstance(widget, cls):
            return widget
        for child in widget.winfo_children():
            found = _find(child, cls)
            if found is not None:
                return found
        return None

    combo = _find(frame, ttk.Combobox)
    assert combo is not None
    assert app.state_.seg_cfg.mode == "auto"

    combo.set("Pathology (no myelin)")
    combo.event_generate("<<ComboboxSelected>>")
    app.update()

    assert app.state_.seg_cfg.mode == "pathological"


def test_gui2_advanced_accordion_preserves_mode_on_apply(app):
    """_build_advanced_accordion's Apply button reconstructs
    SegmentationConfig from only the reflected (non-mode) fields --
    regression test for a real bug caught before commit: without
    explicitly carrying `mode` forward, Apply silently reset it back to
    the dataclass default ("auto")."""
    import dataclasses
    import tkinter as tk
    from tkinter import ttk

    from gui.panels import metrics as metrics_panel

    app.state_.seg_cfg = dataclasses.replace(app.state_.seg_cfg, mode="pathological")

    frame = tk.Frame(app)
    metrics_panel.build(frame, app.state_)

    def _find_button(widget, text):
        if isinstance(widget, ttk.Button) and widget.cget("text") == text:
            return widget
        for child in widget.winfo_children():
            found = _find_button(child, text)
            if found is not None:
                return found
        return None

    apply_btn = _find_button(frame, "Apply")
    assert apply_btn is not None
    apply_btn.invoke()

    assert app.state_.seg_cfg.mode == "pathological"


def test_gui2_maybe_drop_mode_column():
    import pandas as pd

    from gui.panels.export import _maybe_drop_mode_column

    df = pd.DataFrame({"axon_id": [1, 2], "mode": ["normal", "pathological"], "group": ["a", "a"]})

    dropped = _maybe_drop_mode_column(df, include_mode_column=False)
    assert "mode" not in dropped.columns
    assert "group" in dropped.columns

    kept = _maybe_drop_mode_column(df, include_mode_column=True)
    assert "mode" in kept.columns


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

    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    _run_and_wait(app, pairs, "normal", tmp_path)
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

    # Default include_mode_column=False -- mode should be absent.
    import pandas as pd
    written = pd.read_csv(tmp_path / "axons.csv")
    assert "mode" not in written.columns


def test_gui2_results_dashboard_hides_cards_with_no_result(app, tmp_path):
    """Dashboard gap fix: with no result yet, the summary/plot cards
    must be grid_remove()'d entirely (not shown empty, which is what
    caused the reported gap above the summary cards), and must come
    back once a result exists. winfo_manager() (not winfo_ismapped(),
    which is unreliable for a Frame that was never packed/gridded into
    a mapped window -- verified empirically while writing this test) is
    '' for a grid_remove()'d widget and 'grid' once gridded."""
    import tkinter as tk

    from temptation import discovery

    from gui.panels import results as results_panel

    frame = tk.Frame(app)
    results_panel.build(frame, app.state_)
    app.update_idletasks()

    def _find_cards(widget):
        found = []
        if isinstance(widget, tk.Frame) and str(widget.cget("highlightthickness")) not in ("0",):
            found.append(widget)
        for child in widget.winfo_children():
            found.extend(_find_cards(child))
        return found

    cards = _find_cards(frame)
    assert cards, "expected at least the Summary/Group comparison card frames"
    for card in cards:
        assert card.winfo_manager() == ""

    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    _run_and_wait(app, pairs, "normal", tmp_path)
    app.update_idletasks()

    for card in cards:
        assert card.winfo_manager() == "grid"


def test_gui2_review_panel_loads_image_and_inspects_axon(app, tmp_path):
    """Build the review panel directly, select an image (driving its
    on-demand pipeline.analyze_image_legacy call and matplotlib
    overlay render), then simulate an axon click via the panel's own
    hit-testing to confirm the inspector text updates -- exercising the
    on-demand per-image path analyze_dataset deliberately doesn't
    retain (Phase 7 prep's memory tradeoff)."""
    import tkinter as tk

    from temptation import discovery

    from gui.panels import review as review_panel
    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.pairs = pairs
    app.state_.groups["normal"] = entry
    app.state_.pixel_size_um = 0.00524

    frame = tk.Frame(app)
    review_panel.build(frame, app.state_)

    def _find(widget, cls):
        if isinstance(widget, cls):
            return widget
        for child in widget.winfo_children():
            found = _find(child, cls)
            if found is not None:
                return found
        return None

    listbox = _find(frame, tk.Listbox)
    assert listbox is not None
    assert listbox.size() == 4  # 162-165

    # <<ListboxSelect>> only fires from a real mouse click, not from
    # selection_set()/event_generate() outside a live mainloop -- drive
    # the same loader the click handler calls instead (see review.py's
    # parent._load_and_render test hook).
    frame._load_and_render(pairs[0])
    app.update()

    # Inspector text lives in a Label inside a Card widget on the right;
    # search generically rather than depending on the widget tree shape.
    all_text = []

    def _collect_labels(widget):
        if isinstance(widget, tk.Label):
            all_text.append(widget.cget("text"))
        for child in widget.winfo_children():
            _collect_labels(child)

    _collect_labels(frame)
    assert any("axon(s) in" in t for t in all_text)


def test_gui2_review_panel_previous_next_and_list_stays_visible(app):
    """Regression test for the reported bug: the image list appeared to
    disappear once an overlay loaded (an unweighted-column layout bug,
    not a selection bug), and Previous/Next give a layout-independent
    way to change images. Drives real button clicks through a live
    mainloop (not app.update() in a tight loop -- see _run_and_wait's
    docstring; the same segfault risk applies to matplotlib canvas
    creation, confirmed while writing this test)."""
    import tkinter as tk
    from tkinter import ttk

    from temptation import discovery

    from gui.panels import review as review_panel
    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.pairs = pairs
    app.state_.groups["normal"] = entry
    app.state_.pixel_size_um = 0.00524

    top = tk.Toplevel(app)
    frame = tk.Frame(top)
    frame.pack(fill="both", expand=True)
    review_panel.build(frame, app.state_)

    def _find(widget, cls):
        if isinstance(widget, cls):
            return widget
        for child in widget.winfo_children():
            found = _find(child, cls)
            if found is not None:
                return found
        return None

    def _find_button_by_text(widget, text):
        if isinstance(widget, ttk.Button) and widget.cget("text") == text:
            return widget
        for child in widget.winfo_children():
            found = _find_button_by_text(child, text)
            if found is not None:
                return found
        return None

    listbox = _find(frame, tk.Listbox)
    results = {}

    def step0():
        frame._load_and_render(pairs[0])
        app.after(200, step1)

    def step1():
        list_outer = listbox.master.master
        results["list_width_after_load"] = list_outer.winfo_width()
        _find_button_by_text(frame, "Next ▶").invoke()
        app.after(200, step2)

    def step2():
        _find_button_by_text(frame, "Next ▶").invoke()
        app.after(200, step3)

    def step3():
        _find_button_by_text(frame, "◀ Previous").invoke()
        app.after(200, finish)

    def finish():
        texts = []

        def collect(w):
            if isinstance(w, tk.Label):
                texts.append(w.cget("text"))
            for c in w.winfo_children():
                collect(c)

        collect(frame)
        for t in texts:
            if "axon(s) in" in t:
                results["final_inspector_text"] = t
        app.quit()

    app.after(100, step0)
    app.mainloop()

    # minsize configured in review.py's build() is 190px for the list
    # column; a meaningfully larger observed width confirms the canvas
    # isn't crowding it out.
    assert results["list_width_after_load"] >= 190
    # 162 -> Next -> 163 -> Next -> 164 -> Previous -> 163 (14 axons,
    # matches the known-good count from earlier phases' golden checks).
    assert "163" in results["final_inspector_text"]
    assert "14 axon(s)" in results["final_inspector_text"]


def test_gui2_metrics_and_qc_panels_build_with_prior_run_data(app):
    """Building the metrics/QC panels with a populated dataset and a
    prior run's raw axons exercises the mask-scan setup and the QC
    live-preview recompute path, both of which only run once real data
    exists -- a bare construction test wouldn't reach either."""
    import tkinter as tk

    import pandas as pd
    from temptation import discovery
    from temptation.pipeline import DatasetResult

    from gui.panels import metrics as metrics_panel
    from gui.panels import qc as qc_panel
    from gui.state import GroupEntry

    entry = GroupEntry(name="normal")
    pairs, _ = discovery._scan_folder_recursive(NORMAL_DIR)
    entry.pairs = pairs
    app.state_.groups["normal"] = entry
    app.state_.result = DatasetResult(
        df_axons=pd.read_csv(Path(__file__).parent / "golden_v2" / "normal" / "axons.csv")
    )

    metrics_panel.build(tk.Frame(app), app.state_)
    qc_panel.build(tk.Frame(app), app.state_)


def test_gui2_debounce_delays_calibration_commit(app):
    """Calibration's pixel-size field commits state.mark_inputs_changed()
    on a debounced timer, not per keystroke -- regression test for the
    dominant cause of the reported Setup-tab sluggishness (every
    keystroke used to null state.result and fan out to every listener)."""
    import tkinter as tk
    from tkinter import ttk

    from gui.panels import calibration as calibration_panel

    # A real <KeyRelease> (unlike a virtual event such as
    # <<ComboboxSelected>>) needs a mapped, focused widget to dispatch --
    # a bare unpacked tk.Frame(app) doesn't cut it (verified empirically:
    # the event silently never fires and state.pixel_size_um stays None).
    top = tk.Toplevel(app)
    frame = tk.Frame(top)
    frame.pack(fill="both", expand=True)
    calibration_panel.build(frame, app.state_)

    def _find(widget, cls):
        if isinstance(widget, cls):
            return widget
        for child in widget.winfo_children():
            found = _find(child, cls)
            if found is not None:
                return found
        return None

    entry = _find(frame, ttk.Entry)
    assert entry is not None

    # Pre-set a non-default state.step_status so a debounced demotion is
    # observable -- "not_started" is also the AppState.__init__ default,
    # so leaving it untouched wouldn't distinguish "debounce hasn't
    # fired yet" from "there was never anything to fire".
    app.state_.step_status["review"] = "ready"

    results = {}

    def type_value():
        entry.focus_force()
        entry.delete(0, tk.END)
        entry.insert(0, "0.00524")
        entry.event_generate("<KeyRelease>")
        app.after(0, check_immediate)

    def check_immediate():
        # state.pixel_size_um updates immediately (local, cheap)...
        results["pixel_size_um"] = app.state_.pixel_size_um
        # ...but the debounced notify (and step_status demotion) should
        # not have landed yet, right after the keystroke with a 400ms delay.
        results["immediate"] = app.state_.step_status["review"]
        app.after(500, check_settled)

    def check_settled():
        results["settled"] = app.state_.step_status["review"]
        app.quit()

    app.after(50, type_value)
    app.mainloop()

    assert results["pixel_size_um"] == pytest.approx(0.00524)
    assert results["immediate"] == "ready"  # debounce hadn't fired yet
    assert results["settled"] == "not_started"  # fired after the delay
