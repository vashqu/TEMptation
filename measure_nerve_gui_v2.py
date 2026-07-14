#!/usr/bin/env python3

"""
measure_nerve_gui_v2.py
------------------------
Temporary launcher for the in-progress Phase 7 GUI rebuild (gui/app.py,
IMPLEMENTATION_BLUEPRINT.md Sec 8-9). Once gui/ covers the full
workflow (Phase 7d), this replaces measure_nerve_gui.py as the single
launcher; until then both coexist so the old GUI stays usable while the
new one is built incrementally.

Usage:
    python measure_nerve_gui_v2.py
"""

import sys
from pathlib import Path

_this_dir = Path(__file__).resolve().parent
if str(_this_dir) not in sys.path:
    sys.path.insert(0, str(_this_dir))

from gui.app import main  # noqa: E402

if __name__ == "__main__":
    main()
