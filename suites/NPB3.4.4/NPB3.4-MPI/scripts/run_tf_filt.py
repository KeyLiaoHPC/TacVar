#!/usr/bin/env python3
"""
@file run_tf_filt.py
@brief Per (region_id, loc_id): tf = elapsed_tf - ngauge * nspg, then filt.x.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SUITE = Path(__file__).resolve().parents[1]
FILT_SRC = SUITE / "common" / "filt.c"
FILT_HDR = SUITE / "common" / "filt.h"
FILT_EXE = SUITE / "bin" / "filt.x"
DIR_RE = re.compile(r"^([A-Za-z0-9]+)\.([A-Za-z])$")
NGAUGE_MIN = 10
NCDF = 1000


def load_elapsed(dir_path: Path) -> dict[tuple[int, int], list[int]]:
    out: dict[tuple[int, int], list[int]] = defaultdict(list)
    if not dir_path.is_dir():
        return out
    for path in sorted(dir_path.glob("*.csv")):
        if path.name in ("timer_info.csv", "median.csv"):
            continue
        with path.open(newline="") as fp:
            reader = csv.DictReader(fp)
            if not reader.fieldnames or "elapsed_ns" not in reader.fieldnames:
                continue
            for rec in reader:
                try:
                    rid = int(rec["region_id"])
                    loc = int(rec["loc_id"])
                    el = int(float(rec["elapsed_ns"]))
                except (KeyError, ValueError):
                    continue
                out[(rid, loc)].append(el)
    return out


def load_nspg(path: Path) -> float:
    if not path.is_file():
        raise SystemExit(f"missing nspg file: {path}")
    text = path.read_text().strip().splitlines()
    if not text:
        raise SystemExit(f"empty nspg file: {path}")
    nspg = float(text[0].split()[0])
    if nspg <= 0:
        raise SystemExit(f"nspg must be > 0, got {nspg}")
    return nspg


def load_ngauge(path: Path) -> dict[tuple[str, str, int, int], int]:
    table: dict[tuple[str, str, int, int], int] = {}
    if not path.is_file():
        raise SystemExit(f"missing median.csv: {path}")
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        if not reader.fieldnames:
            raise SystemExit(f"empty median.csv: {path}")
        fields = {n.lower(): n for n in reader.fieldnames}
        if "ngauge" not in fields:
            raise SystemExit("median.csv missing ngauge column")
        loc_key = fields.get("loc_id") or fields.get("log_id")
        if loc_key is None:
            raise SystemExit("median.csv missing loc_id")
        for rec in reader:
            try:
                kernel = rec[fields["kernel"]].strip().upper()
                klass = rec[fields["class"]].strip().upper()
                rid = int(rec[fields["region_id"]])
                loc = int(rec[loc_key])
                ngauge = int(float(rec[fields["ngauge"]]))
            except (KeyError, ValueError):
                continue
            table[(kernel, klass, rid, loc)] = ngauge
    return table


def bin_width(arr: np.ndarray) -> int:
    gap = float(np.quantile(arr, 0.5) - np.quantile(arr, 0))
    bw = max(gap / 50.0, 10.0)
    bw = int(bw / 10.0) * 10
    return max(bw, 10)


def ensure_filt_x() -> Path:
    if FILT_EXE.is_file():
        return FILT_EXE
    if not FILT_SRC.is_file() or not FILT_HDR.is_file():
        raise SystemExit(f"missing FilT sources under {SUITE / 'common'}")
    FILT_EXE.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["cc", "-O2", "-std=c99", "-Wall", "-o", str(FILT_EXE), str(FILT_SRC), "-lm"]
    print("building", " ".join(cmd))
    subprocess.check_call(cmd)
    return FILT_EXE


def write_col(path: Path, values) -> None:
    with path.open("w") as fp:
        for v in values:
            fp.write(f"{int(v)}\n")


def write_tf_cdf(path: Path, tf: np.ndarray) -> None:
    qs = np.arange(1, NCDF + 1, dtype=np.float64) / float(NCDF)
    vals = np.quantile(tf, qs)
    with path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["percentile", "tf_ns"])
        for i, v in enumerate(vals, start=1):
            w.writerow([f"{i / 10.0:.1f}", int(np.rint(v))])


def process_pair(
    filt: Path,
    out_dir: Path,
    met: np.ndarray,
    tf: np.ndarray,
    nsamp_sim: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    write_col(out_dir / "met.csv", met)
    write_col(out_dir / "tf.csv", tf)
    write_tf_cdf(out_dir / "tf_cdf.csv", tf)
    bw = bin_width(met)
    cmd = [
        str(filt),
        "-m", "met.csv",
        "-s", "tf.csv",
        "-w", str(bw),
        "-n", str(nsamp_sim),
        "-l", "0.01",
    ]
    subprocess.check_call(cmd, cwd=str(out_dir))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root")
    ap.add_argument("--kernel-class", default="", help="e.g. cg.C (default: all)")
    ap.add_argument("--filt-nsamp", type=int, default=100000)
    args = ap.parse_args()
    root = Path(args.data_root).resolve()
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    names = []
    if args.kernel_class:
        names.append(args.kernel_class)
    else:
        for d in sorted(root.iterdir()):
            if d.is_dir() and DIR_RE.match(d.name):
                names.append(d.name)

    nspg = load_nspg(root / "nspg.txt")
    ngauge_tbl = load_ngauge(root / "median.csv")
    filt = ensure_filt_x()
    n_ok = 0
    n_skip = 0
    for name in names:
        m = DIR_RE.match(name)
        if not m:
            print(f"skip {name}: not Kernel.CLASS", file=sys.stderr)
            n_skip += 1
            continue
        kernel, klass = m.group(1).upper(), m.group(2).upper()
        met_dir = root / name
        tf_dir = root / f"{name}_tf"
        filt_root = root / f"{name}_filt"
        if not met_dir.is_dir() or not tf_dir.is_dir():
            print(f"skip {name}: missing met/tf dirs", file=sys.stderr)
            n_skip += 1
            continue
        met_map = load_elapsed(met_dir)
        tf_map = load_elapsed(tf_dir)
        keys = sorted(set(met_map) & set(tf_map))
        for rid, loc in keys:
            raw_tf = np.asarray(tf_map[(rid, loc)], dtype=np.int64)
            met = np.asarray(met_map[(rid, loc)], dtype=np.int64)
            ngauge = ngauge_tbl.get((kernel, klass, rid, loc), 0)
            if ngauge < NGAUGE_MIN or raw_tf.size == 0 or met.size == 0:
                n_skip += 1
                continue
            tf = raw_tf.astype(np.float64) - float(ngauge) * nspg
            out = filt_root / f"r{rid}_l{loc}"
            try:
                process_pair(filt, out, met, np.rint(tf).astype(np.int64), args.filt_nsamp)
                n_ok += 1
            except subprocess.CalledProcessError as e:
                print(f"filt failed {out}: {e}", file=sys.stderr)
                n_skip += 1
    print(f"run_tf_filt: ok={n_ok} skip={n_skip}")
    if n_ok == 0:
        raise SystemExit("no (region, loc) pairs filtered")


if __name__ == "__main__":
    main()
