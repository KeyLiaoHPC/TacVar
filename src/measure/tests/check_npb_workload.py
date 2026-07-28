#!/usr/bin/env python3
"""
@file check_npb_workload.py
@brief Verify NPB-MPI sources differ from baseline only by TacVar instrumentation.

Strips contiguous TacVar marker blocks:
  Fortran:  !----- TacVar ...   through  !----- end TacVar -----
  C:        /*----- TacVar ...  through  /*----- end TacVar -----*/
then compares SHA-256 of the remaining content to npb_workload_baseline.sha256.
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASELINE = Path(__file__).resolve().parent / "npb_workload_baseline.sha256"

BEGIN_RE = re.compile(r"^\s*(?:!|/\*)-+\s*TacVar\b", re.IGNORECASE)
END_RE = re.compile(r"^\s*(?:!|/\*)-+\s*end\s+TacVar\b", re.IGNORECASE)
LONE_CALL_RE = re.compile(
    r"^\s*(?:call\s+)?tacvar_(?:npb_)?(?:prepare|step_start|step_stop|"
    r"kernel_start|kernel_stop|set_test_tag)\b",
    re.IGNORECASE,
)


def strip_instrumentation(text: str) -> str:
    lines = text.splitlines(keepends=True)
    out: list[str] = []
    skipping = False
    for line in lines:
        if skipping:
            if END_RE.search(line):
                skipping = False
            continue
        if BEGIN_RE.search(line):
            if END_RE.search(line):
                continue
            skipping = True
            continue
        if LONE_CALL_RE.search(line):
            continue
        out.append(line)
    if skipping:
        raise SystemExit("unclosed TacVar instrumentation block")
    # Ensure no residual TacVar API calls remain.
    joined = "".join(out)
    if re.search(r"\btacvar_", joined, re.IGNORECASE):
        raise SystemExit("residual tacvar_ symbols after strip")
    return joined


def load_baseline() -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in BASELINE.read_text().splitlines():
        if not line.strip():
            continue
        digest, path = line.split(None, 1)
        mapping[path] = digest
    return mapping


def main() -> int:
    baseline = load_baseline()
    failed = 0
    for rel, expected in baseline.items():
        path = ROOT / rel
        if not path.exists():
            print(f"MISSING {rel}")
            failed += 1
            continue
        raw_bytes = path.read_bytes()
        raw_hash = hashlib.sha256(raw_bytes).hexdigest()
        if raw_hash == expected:
            print(f"OK unchanged {rel}")
            continue
        try:
            stripped = strip_instrumentation(raw_bytes.decode("utf-8", errors="replace"))
        except SystemExit as exc:
            print(f"FAIL {rel}: {exc}")
            failed += 1
            continue
        stripped_hash = hashlib.sha256(stripped.encode()).hexdigest()
        if stripped_hash != expected:
            print(f"FAIL workload changed beyond TacVar instrumentation: {rel}")
            print(f"  expected {expected}")
            print(f"  got      {stripped_hash}")
            failed += 1
        else:
            print(f"OK instrumented-only {rel}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
