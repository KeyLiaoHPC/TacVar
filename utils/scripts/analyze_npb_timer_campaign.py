#!/usr/bin/env python3
"""Summarize, plot, and write REPORT.md for an NPB-MPI timer campaign."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from npb_timer_campaign.report import analyze_campaign


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results_dir", type=Path)
    p.add_argument("--kernels", default="")
    p.add_argument("--timers", default="")
    p.add_argument("--profiles", default="")
    p.add_argument("--formats", default="png,pdf,eps")
    p.add_argument("--allow-partial", action="store_true")
    p.add_argument("--allow-font-fallback", action="store_true")
    p.add_argument("--full-tail", action="store_true")
    args = p.parse_args(argv)
    kernels = [x for x in args.kernels.split(",") if x] or None
    timers = [x for x in args.timers.split(",") if x] or None
    profiles = [x for x in args.profiles.split(",") if x] or None
    formats = tuple(x.strip() for x in args.formats.split(",") if x.strip())
    out = analyze_campaign(
        args.results_dir,
        formats=formats,
        allow_partial=args.allow_partial,
        allow_font_fallback=args.allow_font_fallback,
        kernels=kernels,
        timers=timers,
        profiles=profiles,
        full_tail=args.full_tail,
    )
    print(out["report"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
