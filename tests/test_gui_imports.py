"""GUI imports must succeed headlessly (no Tk() construction here)."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")


def test_gui_module_imports():
    """measure_nerve_gui.py is now the thin gui/app.py launcher
    (Phase 7d); this only checks the import chain, not construction --
    see tests/test_gui2_shell.py for real Tk construction/behavior."""
    import measure_nerve_gui  # noqa: F401


def test_gui_legacy_module_imports():
    """Retired pre-Phase-7 GUI, kept importable during the transition
    window (blueprint Phase 7 risk note)."""
    import measure_nerve_gui_legacy  # noqa: F401
