"""Shared UI state (IMPLEMENTATION_BLUEPRINT.md Sec 8.11): panels read
and mutate this directly, then call notify() so app.py can refresh the
step rail. Deliberately not a frozen dataclass or an event bus -- a
single-window app with under ten panels and one background worker
thread (marshalled back via Tk's .after()) doesn't need one, and a
plain mutable object keeps each panel's code close to the widget code
it corresponds to."""

import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from temptation.config import QCThresholds, SegmentationConfig

STEP_IDS = ["setup", "review", "export"]

STEP_LABELS = {
    "setup": "① Setup",
    "review": "② Review",
    "export": "③ Export",
}

# Any edit to an input (data/calibration/segmentation/QC, all now living
# inside the single "setup" step) marks a previous run's results stale
# (blueprint Sec 8.1's "results can never be attributed to the wrong
# parameters", adapted after the 7-step -> 3-step reorg: rather than
# nulling `result` outright -- which used to blank the dashboard/export
# checklist mid-keystroke, see AppState.mark_inputs_changed below --
# this only demotes the rail glyph and flips `result_stale`).
DOWNSTREAM_OF_INPUTS = ["review", "export"]


@dataclass
class GroupEntry:
    name: str
    folder: Optional[Path] = None
    pairs: list = field(default_factory=list)
    skipped: list = field(default_factory=list)

    @property
    def n_pairs(self) -> int:
        return len(self.pairs)


class AppState:
    def __init__(self):
        self.groups: dict[str, GroupEntry] = {}
        self.output_dir: Optional[Path] = None

        self.pixel_size_um: Optional[float] = None
        self.pixel_space_only: bool = False

        self.seg_cfg: SegmentationConfig = SegmentationConfig()
        self.qc_thresholds: QCThresholds = QCThresholds()
        self.exclude_qc_failed: bool = False
        self.collect_mito: bool = False
        # Default off: `mode` (the per-image auto-detected segmentation
        # algorithm) is easily mistaken for the `group` label -- see
        # gui/panels/metrics.py's _build_mode_selector -- so it's opt-in
        # in exports rather than a column every user has to explain away.
        self.include_mode_column: bool = False

        self.result = None  # temptation.pipeline.DatasetResult, once a run completes
        # True once an input (data/calibration/segmentation/QC) has been
        # edited since `result` was produced. Panels showing `result`
        # (dashboard, export checklist, QC live preview) keep displaying
        # it -- edits no longer null it out -- but show a staleness note
        # so the user knows a fresh Run would reflect the new parameters.
        self.result_stale: bool = False
        self.run_in_progress: bool = False
        self.cancel_event: threading.Event = threading.Event()

        self.step_status: dict[str, str] = {sid: "not_started" for sid in STEP_IDS}
        self._listeners: list[Callable[[], None]] = []

    def on_change(self, callback: Callable[[], None]) -> None:
        self._listeners.append(callback)

    def notify(self) -> None:
        for cb in self._listeners:
            cb()

    def all_pairs(self) -> list:
        """Every discovered pair across every group, each stamped with
        its group's current label (read fresh here rather than cached on
        the pair dict, since a card's name field can be edited after
        the folder was scanned)."""
        pairs = []
        for entry in self.groups.values():
            for p in entry.pairs:
                p = dict(p)
                p["group"] = entry.name
                pairs.append(p)
        return pairs

    def total_pairs(self) -> int:
        return sum(g.n_pairs for g in self.groups.values())

    def mark_inputs_changed(self) -> None:
        for sid in DOWNSTREAM_OF_INPUTS:
            self.step_status[sid] = "not_started"
        if self.result is not None:
            self.result_stale = True
        self.notify()
