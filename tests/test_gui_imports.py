"""GUI import must succeed headlessly (no Tk() construction here)."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")


def test_gui_module_imports():
    import measure_nerve_gui  # noqa: F401
