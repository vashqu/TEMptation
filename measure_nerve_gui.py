#!/usr/bin/env python3

"""
measure_nerve_gui.py
---------------------
Launches the TEMptation GUI: a step-rail scientific analysis
workstation (IMPLEMENTATION_BLUEPRINT.md Sec 8) covering load data,
calibration, metric selection, QC, run, visual review + results
dashboard, and export -- all backed by the same temptation package the
CLI (measure_nerve.py) uses, never a reimplementation of it (enforced
by tests/test_gui_has_no_metrics.py).

The pre-Phase-7 single-image/batch-folder GUI is retired but still
importable as measure_nerve_gui_legacy.py during this transition
window.

Usage:
    python measure_nerve_gui.py
"""

import sys
from pathlib import Path

_this_dir = Path(__file__).resolve().parent
if str(_this_dir) not in sys.path:
    sys.path.insert(0, str(_this_dir))

from gui.app import main  # noqa: E402

if __name__ == "__main__":
    main()
