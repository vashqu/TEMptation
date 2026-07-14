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

STEP_IDS = ["dataset", "calibration", "metrics", "qc", "run", "review", "export"]

STEP_LABELS = {
    "dataset": "① Load data",
    "calibration": "② Calibrate",
    "metrics": "③ Metrics",
    "qc": "④ QC",
    "run": "⑤ Run",
    "review": "⑥ Review",
    "export": "⑦ Export",
}

# Any edit to an upstream step invalidates a previous run's results
# (blueprint Sec 8.1: "Any change to ①-④ demotes ⑤-⑦ to not started,
# and greys the Export button, so results can never be attributed to
# the wrong parameters").
DOWNSTREAM_OF_INPUTS = ["run", "review", "export"]


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

        self.result = None  # temptation.pipeline.DatasetResult, once a run completes
        # Raw per-axon metrics from the last completed run, kept around
        # independent of `result`/mark_inputs_changed() so the QC panel's
        # live preview (blueprint Sec 8.5: "recomputed on threshold edit
        # ... no re-segmentation") still has something to recompute
        # against even after a threshold edit demotes the Run step.
        self.last_raw_axons = None
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
        self.result = None
        self.notify()
