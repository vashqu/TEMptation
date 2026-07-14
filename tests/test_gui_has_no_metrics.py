"""Enforces "no metric logic in the GUI" (IMPLEMENTATION_BLUEPRINT.md
Sec 4.4). Parses GUI source with ast rather than importing it, so this
test runs even without a display. Covers both the legacy
measure_nerve_gui.py and the in-progress gui/ package (Phase 7,
blueprint Sec 9's Phase 7 checklist: "test_gui_has_no_metrics.py
extended over gui/") -- every .py file under either surface must stay
formula-free for the guard to mean anything as the rebuild proceeds."""

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
GUI_FILES = [REPO_ROOT / "measure_nerve_gui.py"] + sorted(
    (REPO_ROOT / "gui").rglob("*.py")
)

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
    # qc.compute_qc_flags is itself the backend flag engine -- the QC
    # panel's live preview (blueprint Sec 8.5) calls it directly so the
    # preview can never drift from what a real Run would compute; that
    # is reuse, not the reimplementation this guard exists to catch.
    "qc",
    # summaries.summarize_groups is the same backend function cli.py
    # calls for --write-group-csv; the export panel calls it to write
    # group_metrics.csv and to estimate its row count for the checklist.
    "summaries",
}


def _called_name(call_node: ast.Call):
    func = call_node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def test_gui_has_no_banned_calls():
    offenders = []
    for path in GUI_FILES:
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                name = _called_name(node)
                if name in BANNED_CALLS:
                    offenders.append((path.relative_to(REPO_ROOT), name, node.lineno))
    assert not offenders, (
        f"GUI code calls segmentation/measurement functions directly "
        f"(should call temptation.pipeline instead): {offenders}"
    )


def test_gui_temptation_imports_are_restricted():
    offenders = []
    for path in GUI_FILES:
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if module == "temptation":
                    # from temptation import X[, Y...]
                    for alias in node.names:
                        if alias.name not in ALLOWED_TEMPTATION_SUBMODULES:
                            offenders.append((path.relative_to(REPO_ROOT), module, alias.name, node.lineno))
                elif module.startswith("temptation."):
                    submodule = module.split(".", 1)[1]
                    if submodule not in ALLOWED_TEMPTATION_SUBMODULES:
                        offenders.append((path.relative_to(REPO_ROOT), module, None, node.lineno))
    assert not offenders, (
        f"GUI code imports a non-approved temptation submodule "
        f"(only {sorted(ALLOWED_TEMPTATION_SUBMODULES)} are allowed): {offenders}"
    )


def test_gui_has_no_pi_arithmetic():
    """Heuristic: a BinOp involving `np.pi` or a bare `pi` name is the
    signature of a hand-rolled circularity/area formula."""
    offenders = []
    for path in GUI_FILES:
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.BinOp):
                for side in (node.left, node.right):
                    if isinstance(side, ast.Attribute) and side.attr == "pi":
                        offenders.append((path.relative_to(REPO_ROOT), node.lineno))
                    if isinstance(side, ast.Name) and side.id == "pi":
                        offenders.append((path.relative_to(REPO_ROOT), node.lineno))
    assert not offenders, f"GUI code contains pi-based arithmetic at: {offenders}"
