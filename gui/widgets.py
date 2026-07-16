"""Design tokens and small reusable widgets (cards, scrollable areas, tooltips, etc.) for the TEMptation GUI"""

import tkinter as tk
from tkinter import ttk

BG          = "#F5F7FA"
CARD_BG     = "white"
HDR_FG      = "#333333"
LBL_FG      = "#666666"
BORDER_CLR  = "#CCCCCC"
BTN_PRIMARY = "#007BFF"
BTN_HOVER   = "#0056b3"
BTN_FG      = "white"
STATUS_BG   = "#EFEFEF"
EMPTY_FG    = "#888888"
ROW_ODD     = "white"
ROW_EVEN    = "#F9F9F9"

WARN_FG    = "#B7791F"
ERROR_FG   = "#C53030"
SUCCESS_FG = "#2F855A"

CARD_PAD = 14

# One glyph + color per step-rail status (blueprint Sec 8.1).
STATUS_GLYPHS = {
    "not_started": ("○", LBL_FG),
    "ready":       ("●", BTN_PRIMARY),
    "complete":    ("✓", SUCCESS_FG),
    "warning":     ("⚠", WARN_FG),
}


class Tooltip:
    """Show a small floating tooltip when the mouse hovers over a widget."""

    def __init__(self, widget: tk.Widget, text: str):
        self._widget = widget
        self._text = text
        self._tip: tk.Toplevel | None = None
        widget.bind("<Enter>", self._show, add="+")
        widget.bind("<Leave>", self._hide, add="+")

    def _show(self, _event=None):
        if self._tip or not self._text:
            return
        x = self._widget.winfo_rootx() + 20
        y = self._widget.winfo_rooty() + self._widget.winfo_height() + 4
        self._tip = tw = tk.Toplevel(self._widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        tk.Label(
            tw,
            text=self._text,
            justify="left",
            background="#ffffc0",
            foreground="#222222",
            relief="solid",
            borderwidth=1,
            font=("TkDefaultFont", 9),
            wraplength=340,
            padx=6,
            pady=5,
        ).pack()

    def _hide(self, _event=None):
        if self._tip:
            self._tip.destroy()
            self._tip = None


def make_card(parent, title: str, **grid_kw) -> tuple[tk.Frame, tk.Frame]:
    """White card with a bold title. Returns (outer_frame, inner_frame)."""
    outer = tk.Frame(parent, bg=CARD_BG,
                      highlightthickness=1, highlightbackground=BORDER_CLR)
    if grid_kw:
        outer.grid(**grid_kw)
    if title:
        tk.Label(outer, text=title, bg=CARD_BG, fg=HDR_FG,
                  font=("TkDefaultFont", 11, "bold")).pack(
            anchor="w", padx=CARD_PAD, pady=(CARD_PAD, 4))
    inner = tk.Frame(outer, bg=CARD_BG)
    inner.pack(fill="both", expand=True, padx=CARD_PAD, pady=(0, CARD_PAD))
    return outer, inner


def make_scrollable(parent, bg: str = None) -> tuple[tk.Frame, tk.Frame]:
    """Wraps content in a vertically scrollable area (Canvas + Scrollbar
    + an inner Frame that holds the real content). Build into the
    returned inner frame exactly as you would a plain Frame; grid/pack
    the returned outer frame into `parent`.

    Used for the Setup step's sub-tabs: their combined content (dataset
    cards, ~15 segmentation parameters, ~7 QC thresholds, seven metric
    cards) can exceed the window height with no way to reach the rest,
    e.g. the "Image-level summaries" card being cut off at the bottom
    with nothing below it visible or reachable."""
    bg = bg or BG
    outer = tk.Frame(parent, bg=bg)
    outer.rowconfigure(0, weight=1)
    outer.columnconfigure(0, weight=1)

    canvas = tk.Canvas(outer, bg=bg, highlightthickness=0)
    canvas.grid(row=0, column=0, sticky="nsew")
    scrollbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
    scrollbar.grid(row=0, column=1, sticky="ns")
    canvas.configure(yscrollcommand=scrollbar.set)

    inner = tk.Frame(canvas, bg=bg)
    inner_window = canvas.create_window((0, 0), window=inner, anchor="nw")

    def _on_inner_configure(_event=None):
        canvas.configure(scrollregion=canvas.bbox("all"))

    def _on_canvas_configure(event):
        # inner's width tracks the canvas's -- only its height should
        # ever exceed the viewport, never its width (no horizontal
        # scrollbar is offered).
        canvas.itemconfigure(inner_window, width=event.width)

    inner.bind("<Configure>", _on_inner_configure)
    canvas.bind("<Configure>", _on_canvas_configure)

    def _on_mousewheel(event):
        if event.num == 5 or event.delta < 0:
            canvas.yview_scroll(1, "units")
        elif event.num == 4 or event.delta > 0:
            canvas.yview_scroll(-1, "units")

    def _bind_mousewheel(_event=None):
        # bind_all is scoped to hover (Enter/Leave), not permanent, so
        # scrolling elsewhere in the app is never hijacked by whichever
        # scrollable tab was last visited.
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        canvas.bind_all("<Button-4>", _on_mousewheel)
        canvas.bind_all("<Button-5>", _on_mousewheel)

    def _unbind_mousewheel(_event=None):
        canvas.unbind_all("<MouseWheel>")
        canvas.unbind_all("<Button-4>")
        canvas.unbind_all("<Button-5>")

    canvas.bind("<Enter>", _bind_mousewheel)
    canvas.bind("<Leave>", _unbind_mousewheel)

    return outer, inner


def debounce(widget: tk.Misc, delay_ms: int, callback):
    """Returns a wrapper that delays `callback` by delay_ms, canceling
    any pending call from a previous invocation. Bind the wrapper to a
    per-keystroke event (KeyRelease, a StringVar trace) instead of
    `callback` directly so rapid typing collapses into one call after a
    pause, rather than firing state.notify() -- and every listener it
    cascades to -- on every keystroke (the dominant cause of the GUI's
    Setup-tab sluggishness before this existed)."""
    pending = {"id": None}

    def _debounced(*args, **kwargs):
        if pending["id"] is not None:
            widget.after_cancel(pending["id"])
        pending["id"] = widget.after(delay_ms, lambda: callback(*args, **kwargs))

    return _debounced


def apply_base_styles(root: tk.Tk) -> ttk.Style:
    s = ttk.Style(root)
    try:
        s.theme_use("clam")
    except Exception:
        pass

    s.configure("TFrame",      background=BG)
    s.configure("Card.TFrame", background=CARD_BG)

    s.configure("TLabel",      background=BG,      foreground=LBL_FG)
    s.configure("Card.TLabel", background=CARD_BG, foreground=LBL_FG)

    s.configure("TEntry",
                fieldbackground="white",
                bordercolor=BORDER_CLR,
                lightcolor=BORDER_CLR,
                darkcolor=BORDER_CLR,
                insertcolor=HDR_FG)

    s.configure("Ghost.TButton",
                background=CARD_BG,
                foreground="#444444",
                borderwidth=1,
                relief="solid",
                padding=(8, 4))
    s.map("Ghost.TButton",
          background=[("active", "#F0F0F0"), ("disabled", CARD_BG)],
          foreground=[("disabled", "#AAAAAA")])

    s.configure("TCombobox",
                fieldbackground="white",
                bordercolor=BORDER_CLR)

    s.configure("Treeview",
                background=ROW_ODD,
                fieldbackground=ROW_ODD,
                rowheight=22,
                foreground=HDR_FG)
    s.configure("Treeview.Heading",
                background="#EAECEF",
                foreground=HDR_FG,
                relief="flat",
                font=("TkDefaultFont", 9, "bold"))
    s.map("Treeview",
          background=[("selected", BTN_PRIMARY)],
          foreground=[("selected", "white")])

    s.configure("Card.TCheckbutton", background=CARD_BG, foreground=LBL_FG)

    s.configure("TProgressbar", background=BTN_PRIMARY, troughcolor="#DDE1E7")

    return s
