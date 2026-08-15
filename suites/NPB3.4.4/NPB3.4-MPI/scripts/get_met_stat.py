#!/usr/bin/env python3
"""
@file get_met_stat.py
@brief Per-site median/p95 of elapsed_ns and ngauge = round(median/nspg).
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


def load_timer_event(data_root: Path) -> tuple[str, str]:
    timer = "unknown"
    event = "none"
    info = data_root / "timer_info.csv"
    if info.is_file():
        with info.open(newline="") as fp:
            reader = csv.DictReader(fp)
            for rec in reader:
                if rec.get("timer"):
                    timer = rec["timer"].strip()
                    break
    for _k, _c, path in iter_rank_csvs(data_root):
        with path.open(newline="") as fp:
            reader = csv.DictReader(fp)
            if reader.fieldnames and "counter_backend" in reader.fieldnames:
                for rec in reader:
                    ev = (rec.get("counter_backend") or "").strip()
                    if ev:
                        event = ev
                        return timer, event
            break
    return timer, event


def load_nspg_file(path: Path) -> float | None:
    if not path.is_file():
        return None
    text = path.read_text().split()
    if not text:
        return None
    try:
        v = float(text[0])
    except ValueError:
        return None
    return v if v > 0.0 else None


def load_nspg(data_root: Path) -> float | None:
    return load_nspg_file(data_root / "nspg.txt")


def read_elapsed(path: Path) -> list[tuple[int, int, int, int]]:
    """(region_id, loc_id, elapsed_ns, rank)"""
    rows: list[tuple[int, int, int, int]] = []
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        if not reader.fieldnames or "elapsed_ns" not in reader.fieldnames:
            return rows
        for rec in reader:
            try:
                rid = int(rec["region_id"])
                loc = int(rec["loc_id"])
                el = int(float(rec["elapsed_ns"]))
                rank = int(float(rec.get("rank", "0") or 0))
            except (KeyError, ValueError):
                continue
            rows.append((rid, loc, el, rank))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", help="data_YYYYMMDDTHHmmss or combo/met directory")
    ap.add_argument("--nspg", default="", help="path to nspg.txt (default: data_root/nspg.txt)")
    args = ap.parse_args()
    root = Path(args.data_root).resolve()
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    timer, event = load_timer_event(root)
    nspg = load_nspg_file(Path(args.nspg)) if args.nspg else load_nspg(root)
    buckets: dict[tuple[str, str, int, int], list[int]] = defaultdict(list)
    ranks: dict[tuple[str, str, int, int], set[int]] = defaultdict(set)
    nfiles = 0
    for kernel, klass, path in iter_rank_csvs(root):
        nfiles += 1
        for rid, loc, el, rank in read_elapsed(path):
            key = (kernel, klass, rid, loc)
            buckets[key].append(el)
            ranks[key].add(rank)

    if not buckets:
        raise SystemExit(f"no measurement CSVs under {root} (looked for Kernel.CLASS/)")

    met_path = root / "met_stat.csv"
    med_path = root / "median.csv"
    with met_path.open("w", newline="") as fp_m, med_path.open("w", newline="") as fp_d:
        wm = csv.writer(fp_m, lineterminator="\n")
        wd = csv.writer(fp_d, lineterminator="\n")
        wm.writerow(
            [
                "kernel",
                "class",
                "timer_reader",
                "event_reader",
                "region_id",
                "loc_id",
                "nmet_total",
                "nrank",
                "p50",
                "p95",
            ]
        )
        wd.writerow(["kernel", "class", "region_id", "loc_id", "median", "ngauge"])
        for key in sorted(buckets):
            arr = np.asarray(buckets[key], dtype=np.int64)
            p50 = int(np.median(arr))
            p95 = int(np.percentile(arr, 95))
            nrank = len(ranks[key])
            ngauge = 0
            if nspg is not None and nspg > 0.0:
                ngauge = int(p50 / nspg + 0.5)
            wm.writerow(
                [
                    key[0],
                    key[1],
                    timer,
                    event,
                    key[2],
                    key[3],
                    len(arr),
                    nrank,
                    p50,
                    p95,
                ]
            )
            wd.writerow([key[0], key[1], key[2], key[3], p50, ngauge])
    print(f"wrote {met_path} and {med_path} ({len(buckets)} sites from {nfiles} files)")


if __name__ == "__main__":
    main()
