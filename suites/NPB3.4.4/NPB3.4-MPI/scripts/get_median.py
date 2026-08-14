#!/usr/bin/env python3
"""
@file get_median.py
@brief Median elapsed_ns per (kernel, class, region_id, loc_id) across all ranks.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

DIR_RE = re.compile(r"^([A-Za-z0-9]+)\.([A-Za-z])$")


def iter_rank_csvs(data_root: Path):
    for d in sorted(data_root.iterdir()):
        if not d.is_dir():
            continue
        m = DIR_RE.match(d.name)
        if not m:
            continue
        kernel = m.group(1).upper()
        klass = m.group(2).upper()
        for csv_path in sorted(d.glob("*.csv")):
            if csv_path.name == "timer_info.csv":
                continue
            yield kernel, klass, csv_path


def read_elapsed(path: Path) -> list[tuple[int, int, int]]:
    rows: list[tuple[int, int, int]] = []
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        if not reader.fieldnames or "elapsed_ns" not in reader.fieldnames:
            return rows
        for rec in reader:
            try:
                rid = int(rec["region_id"])
                loc = int(rec["loc_id"])
                el = int(float(rec["elapsed_ns"]))
            except (KeyError, ValueError):
                continue
            rows.append((rid, loc, el))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", help="data_YYYYMMDDTHHmmss directory")
    args = ap.parse_args()
    root = Path(args.data_root).resolve()
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    buckets: dict[tuple[str, str, int, int], list[int]] = defaultdict(list)
    nfiles = 0
    for kernel, klass, path in iter_rank_csvs(root):
        nfiles += 1
        for rid, loc, el in read_elapsed(path):
            buckets[(kernel, klass, rid, loc)].append(el)

    if not buckets:
        raise SystemExit(f"no measurement CSVs under {root} (looked for Kernel.CLASS/)")

    out = root / "median.csv"
    with out.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["kernel", "class", "region_id", "loc_id", "median"])
        for key in sorted(buckets):
            med = int(np.median(np.asarray(buckets[key], dtype=np.int64)))
            w.writerow([key[0], key[1], key[2], key[3], med])
    print(f"wrote {out} ({len(buckets)} sites from {nfiles} files)")


if __name__ == "__main__":
    main()
