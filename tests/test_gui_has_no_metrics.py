"""Enforces "no metric logic in the GUI" (IMPLEMENTATION_BLUEPRINT.md
Sec 4.4). Parses measure_nerve_gui.py with ast rather than importing it,
so this test runs even without a display."""

import ast
from pathlib import Path

GUI_PATH = Path(__file__).parent.parent / "measure_nerve_gui.py"

# Functions that perform actual segmentation/measurement work. If any of
# these show up in the GUI source, formula/logic has leaked out of the
# temptation package.
BANNED_CALLS = {
    "regionprops", "label", "watershed",
    "binary_dilation", "binary_erosion", "binary_opening", "binary_closing",
    "binary_fill_holes", "remove_small_objects",
    "distance_transform_edt", "find_contours",
}

ALLOWED_TEMPTATION_SUBMODULES = {
    "pipeline", "plotting", "export", "config", "dataio", "discovery", "compat",
}


def _called_name(call_node: ast.Call):
    func = call_node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def test_gui_has_no_banned_calls():
    tree = ast.parse(GUI_PATH.read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = _called_name(node)
            if name in BANNED_CALLS:
                offenders.append((name, node.lineno))
    assert not offenders, (
        f"measure_nerve_gui.py calls segmentation/measurement functions "
        f"directly (should call temptation.pipeline instead): {offenders}"
    )


def test_gui_temptation_imports_are_restricted():
    tree = ast.parse(GUI_PATH.read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if module == "temptation":
                # from temptation import X[, Y...]
                for alias in node.names:
                    if alias.name not in ALLOWED_TEMPTATION_SUBMODULES:
                        offenders.append((module, alias.name, node.lineno))
            elif module.startswith("temptation."):
                submodule = module.split(".", 1)[1]
                if submodule not in ALLOWED_TEMPTATION_SUBMODULES:
                    offenders.append((module, None, node.lineno))
    assert not offenders, (
        f"measure_nerve_gui.py imports a non-approved temptation submodule "
        f"(only {sorted(ALLOWED_TEMPTATION_SUBMODULES)} are allowed): {offenders}"
    )


def test_gui_has_no_pi_arithmetic():
    """Heuristic: a BinOp involving `np.pi` or a bare `pi` name is the
    signature of a hand-rolled circularity/area formula."""
    tree = ast.parse(GUI_PATH.read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.BinOp):
            for side in (node.left, node.right):
                if isinstance(side, ast.Attribute) and side.attr == "pi":
                    offenders.append(node.lineno)
                if isinstance(side, ast.Name) and side.id == "pi":
                    offenders.append(node.lineno)
    assert not offenders, f"measure_nerve_gui.py contains pi-based arithmetic at lines: {offenders}"
