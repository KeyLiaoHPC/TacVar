#!/usr/bin/env python3
"""Compatibility wrapper: summarize a Class C (or any) NPB-MPI timer result tree.

Delegates to utils/scripts/analyze_npb_timer_campaign.py. Legacy c920bn1/Class C
defaults apply only when pointing at that historical result directory — they are
not implicit campaign defaults for new runs.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS = _ROOT / "utils" / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from npb_timer_campaign.report import analyze_campaign  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "results_dir",
        nargs="?",
        default="",
        help="Results directory (default: newest results_c920bn1_classC_*)",
    )
    p.add_argument("--allow-partial", action="store_true", default=True)
    p.add_argument("--allow-font-fallback", action="store_true")
    p.add_argument("--formats", default="png,pdf,eps")
    args = p.parse_args(argv)

    results = Path(args.results_dir) if args.results_dir else None
    if results is None:
        mpi = Path(__file__).resolve().parent / "NPB3.4-MPI"
        cands = sorted(mpi.glob("results_c920bn1_classC_*"), key=lambda p: p.stat().st_mtime)
        if not cands:
            print("No results_c920bn1_classC_* found; pass results_dir explicitly.", file=sys.stderr)
            return 2
        results = cands[-1]
        print(f"[wrapper] using {results}", file=sys.stderr)

    formats = tuple(x.strip() for x in args.formats.split(",") if x.strip())
    out = analyze_campaign(
        results,
        formats=formats,
        allow_partial=args.allow_partial,
        allow_font_fallback=args.allow_font_fallback,
    )
    print(out["report"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
