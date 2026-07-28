#!/usr/bin/env python3
"""Plot normalized per-step histograms for an NPB-MPI timer campaign."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from npb_timer_campaign.plots import plot_step_histograms


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results_dir", type=Path)
    p.add_argument("--kernels", default="", help="Comma-separated filter")
    p.add_argument("--timers", default="", help="Comma-separated filter")
    p.add_argument("--profiles", default="", help="Comma-separated filter")
    p.add_argument("--formats", default="png,pdf,eps")
    p.add_argument("--allow-font-fallback", action="store_true")
    p.add_argument("--nbins", type=int, default=24)
    args = p.parse_args(argv)
    kernels = [x for x in args.kernels.split(",") if x] or None
    timers = [x for x in args.timers.split(",") if x] or None
    profiles = [x for x in args.profiles.split(",") if x] or None
    formats = tuple(x.strip() for x in args.formats.split(",") if x.strip())
    out = plot_step_histograms(
        args.results_dir,
        formats=formats,
        kernels=kernels,
        timers=timers,
        profiles=profiles,
        allow_font_fallback=args.allow_font_fallback,
        nbins=args.nbins,
    )
    print(out.get("json", ""))
    for path in out.get("paths") or []:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
