"""Academic matplotlib style: Times New Roman required unless fallback allowed."""
from __future__ import annotations

from typing import Sequence

import matplotlib as mpl
from matplotlib import font_manager


_PREFERRED = "Times New Roman"
_FALLBACKS = ("Times", "Nimbus Roman", "Liberation Serif", "DejaVu Serif")


def times_new_roman_available() -> bool:
    names = {f.name for f in font_manager.fontManager.ttflist}
    return _PREFERRED in names


def resolve_serif_font(allow_font_fallback: bool = False) -> str:
    names = {f.name for f in font_manager.fontManager.ttflist}
    if _PREFERRED in names:
        return _PREFERRED
    if not allow_font_fallback:
        raise RuntimeError(
            "Times New Roman is required for publication figures but was not found "
            "by matplotlib.font_manager. Install a Times New Roman / msttcorefonts "
            "package, or pass allow_font_fallback=True / --allow-font-fallback "
            "(recorded in REPORT.md)."
        )
    for cand in _FALLBACKS:
        if cand in names:
            return cand
    return "serif"


def apply_academic_style(
    *,
    allow_font_fallback: bool = False,
    font_size: int = 9,
) -> dict[str, str]:
    family = resolve_serif_font(allow_font_fallback=allow_font_fallback)
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [family],
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size + 1,
            "xtick.labelsize": font_size - 1,
            "ytick.labelsize": font_size - 1,
            "legend.fontsize": font_size - 1,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.grid": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )
    return {
        "requested_font": _PREFERRED,
        "resolved_font": family,
        "font_fallback": str(family != _PREFERRED).lower(),
    }


# Grayscale-safe cycling styles for overlaying timers
LINE_STYLES: Sequence[str] = ("-", "--", "-.", ":", (0, (3, 1, 1, 1)))
MARKERS: Sequence[str] = ("o", "s", "^", "D", "v", "x", "+", "P")
