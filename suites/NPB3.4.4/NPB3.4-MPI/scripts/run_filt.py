#!/usr/bin/env python3
"""
@file run_filt.py
@brief Combine tfs/tfe gauge samples into tf, extract met, run filt.x per site.
"""
from __future__ import annotations

import argparse
import csv
import random
import re
import subprocess
import sys
from pathlib import Path

DIR_RE = re.compile(r"^([A-Za-z0-9]+)\.([A-Za-z])$")
TF_RE = re.compile(r"^([A-Za-z0-9]+)\.([A-Za-z])_(tfs|tfe)$")
SITE_RE = re.compile(r"^([A-Za-z0-9]+)\.([A-Za-z])_r(\d+)_l(\d+)$")


def is_new_tree(root: Path) -> bool:
    return (root / "met").is_dir() and (root / "tf").is_dir()


def load_nspg(root: Path) -> float:
    candidates = [
        root / "nspg.txt",
        root / "met" / "nspg.txt",
        root.parent / "nspg.txt",
    ]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        raise SystemExit(f"missing nspg.txt under {root} or {root.parent}")
    v = float(path.read_text().split()[0])
    if v <= 0.0:
        raise SystemExit(f"bad nspg in {path}")
    return v


def load_ngauge_table(root: Path) -> dict[tuple[str, str, int, int], tuple[int, int]]:
    """(kernel, class, rid, loc) -> (median, ngauge)."""
    out: dict[tuple[str, str, int, int], tuple[int, int]] = {}
    path = root / "median.csv"
    if not path.is_file():
        path = root / "met" / "median.csv"
    if not path.is_file():
        raise SystemExit(f"missing median.csv under {root} or {root}/met")
    with path.open(newline="") as fp:
        for rec in csv.DictReader(fp):
            try:
                k = rec["kernel"].upper()
                c = rec["class"].upper()
                rid = int(rec["region_id"])
                loc = int(rec["loc_id"])
                med = int(float(rec["median"]))
                ng = int(float(rec.get("ngauge") or 0))
            except (KeyError, ValueError):
                continue
            out[(k, c, rid, loc)] = (med, ng)
    return out


def read_elapsed(path: Path, rid: int | None = None, loc: int | None = None) -> list[int]:
    vals: list[int] = []
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        if not reader.fieldnames or "elapsed_ns" not in reader.fieldnames:
            return vals
        for rec in reader:
            try:
                r = int(rec["region_id"])
                l = int(rec["loc_id"])
                el = int(float(rec["elapsed_ns"]))
            except (KeyError, ValueError):
                continue
            if rec.get("migrated") in ("1", "1.0"):
                continue
            if rid is not None and r != rid:
                continue
            if loc is not None and l != loc:
                continue
            vals.append(el)
    return vals


def collect_site_elapsed(kernel_dir: Path, rid: int, loc: int) -> list[int]:
    out: list[int] = []
    for csv_path in sorted(kernel_dir.glob("*.csv")):
        if csv_path.name == "timer_info.csv":
            continue
        out.extend(read_elapsed(csv_path, rid, loc))
    return out


def find_kernel_dirs(root: Path) -> dict[tuple[str, str], Path]:
    found: dict[tuple[str, str], Path] = {}
    for d in root.iterdir():
        if not d.is_dir():
            continue
        m = DIR_RE.match(d.name)
        if m:
            found[(m.group(1).upper(), m.group(2).upper())] = d
    return found


def find_tf_dirs(root: Path) -> dict[tuple[str, str, str], Path]:
    found: dict[tuple[str, str, str], Path] = {}
    for d in root.iterdir():
        if not d.is_dir():
            continue
        m = TF_RE.match(d.name)
        if m:
            found[(m.group(1).upper(), m.group(2).upper(), m.group(3))] = d
    return found


def write_col(path: Path, vals: list[int]) -> None:
    with path.open("w") as fp:
        for v in vals:
            fp.write(f"{int(v)}\n")


def combine_tf(tfs: list[int], tfe: list[int], nspg: float, ngauge: int) -> list[int]:
    c = nspg * float(ngauge)
    n = min(len(tfs), len(tfe))
    if n < 1:
        return []
    tfe_sh = tfe[:]
    random.shuffle(tfe_sh)
    out = []
    for i in range(n):
        out.append(int(round((tfs[i] - c) + (tfe_sh[i] - c))))
    return out


def find_filt_x(root: Path, explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit)
        if p.is_file():
            return p
        raise SystemExit(f"filt.x not found: {p}")
    walk = root
    for _ in range(6):
        cand = walk / "bin" / "filt.x"
        if cand.is_file():
            return cand
        if walk.parent == walk:
            break
        walk = walk.parent
    candidates = [
        root.parent / "bin" / "filt.x",
        Path.cwd() / "bin" / "filt.x",
        root / "filt.x",
    ]
    for p in candidates:
        if p.is_file():
            return p
    raise SystemExit("filt.x not found (build with make filt, or pass --filt)")


def resolve_ngauge(med: int, ng: int, nspg: float) -> int:
    if ng >= 1:
        return ng
    if med > 0 and nspg > 0.0:
        return int(med / nspg + 0.5)
    return 0


def collect_sites(
    k: str,
    c: str,
    tfs_dir: Path,
    ngauge_tab: dict[tuple[str, str, int, int], tuple[int, int]],
    nspg: float,
    rid_f: int | None,
    lid_f: int | None,
) -> list[tuple[int, int, int]]:
    """Return (rid, loc, ngauge) for this kernel, optionally filtered."""
    cand: dict[tuple[int, int], int] = {}
    for (kk, cc, rid, loc), (med, ng) in ngauge_tab.items():
        if kk != k or cc != c:
            continue
        cand[(rid, loc)] = resolve_ngauge(med, ng, nspg)
    for csv_path in tfs_dir.glob("*.csv"):
        if csv_path.name == "timer_info.csv":
            continue
        with csv_path.open(newline="") as fp:
            reader = csv.DictReader(fp)
            for rec in reader:
                try:
                    rid = int(rec["region_id"])
                    loc = int(rec["loc_id"])
                except (KeyError, ValueError):
                    continue
                if (rid, loc) in cand:
                    continue
                med, ng = ngauge_tab.get((k, c, rid, loc), (0, 0))
                cand[(rid, loc)] = resolve_ngauge(med, ng, nspg)
    out = []
    for (rid, loc), ng in sorted(cand.items()):
        if rid_f is not None and rid != rid_f:
            continue
        if lid_f is not None and loc != lid_f:
            continue
        out.append((rid, loc, ng))
    return out


def collect_new_tree_jobs(
    root: Path,
    ngauge_tab: dict[tuple[str, str, int, int], tuple[int, int]],
    nspg: float,
    rid_f: int | None,
    lid_f: int | None,
    kernel_f: str | None = None,
) -> list[tuple[str, str, Path, Path, Path, Path, int, int, int]]:
    """(k, c, met_dir, tfs_dir, tfe_dir, out_dir, rid, loc, ngauge)."""
    jobs = []
    tf_root = root / "tf"
    met_root = root / "met"
    filt_root = root / "filter"
    kern = kernel_f.upper() if kernel_f else None
    if not tf_root.is_dir():
        return jobs
    for site_dir in sorted(tf_root.iterdir()):
        if not site_dir.is_dir():
            continue
        m = SITE_RE.match(site_dir.name)
        if not m:
            continue
        k = m.group(1).upper()
        c = m.group(2).upper()
        rid = int(m.group(3))
        loc = int(m.group(4))
        if kern is not None and k != kern:
            continue
        if rid_f is not None and rid != rid_f:
            continue
        if lid_f is not None and loc != lid_f:
            continue
        tfs_dir = site_dir / "tfs"
        tfe_dir = site_dir / "tfe"
        met_dir = met_root / f"{m.group(1)}.{m.group(2)}"
        if not met_dir.is_dir():
            met_dir = met_root / f"{k}.{c}"
        if not (tfs_dir.is_dir() and tfe_dir.is_dir() and met_dir.is_dir()):
            continue
        med, ng = ngauge_tab.get((k, c, rid, loc), (0, 0))
        jobs.append(
            (
                k,
                c,
                met_dir,
                tfs_dir,
                tfe_dir,
                filt_root / site_dir.name,
                rid,
                loc,
                resolve_ngauge(med, ng, nspg),
            )
        )
    return jobs


def run_one_site(
    k: str,
    c: str,
    met_dir: Path,
    tfs_dir: Path,
    tfe_dir: Path,
    out_dir: Path,
    rid: int,
    loc: int,
    ng: int,
    nspg: float,
    filt_x: Path,
    width: int,
    nsamp: int,
    plow: float,
) -> bool:
    if (out_dir / "wd.out").is_file():
        print(f"skip {k}.{c} r{rid} l{loc}: wd.out exists", file=sys.stderr)
        return True
    if ng < 1:
        print(f"skip {k}.{c} r{rid} l{loc}: ngauge < 1", file=sys.stderr)
        return False
    tfs = collect_site_elapsed(tfs_dir, rid, loc)
    tfe = collect_site_elapsed(tfe_dir, rid, loc)
    met = collect_site_elapsed(met_dir, rid, loc)
    tf = combine_tf(tfs, tfe, nspg, ng)
    if not met or not tf:
        print(
            f"skip {k}.{c} r{rid} l{loc}: met={len(met)} tf={len(tf)} "
            f"tfs={len(tfs)} tfe={len(tfe)}",
            file=sys.stderr,
        )
        return False
    out_dir.mkdir(parents=True, exist_ok=True)
    write_col(out_dir / "met.csv", met)
    write_col(out_dir / "tf.csv", tf)
    cmd = [
        str(filt_x),
        "-m",
        "met.csv",
        "-s",
        "tf.csv",
        "-w",
        str(width),
        "-n",
        str(nsamp),
        "-l",
        str(plow),
    ]
    print(f"filt {out_dir} (nmet={len(met)} ntf={len(tf)} ngauge={ng})")
    subprocess.check_call(cmd, cwd=out_dir)
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", help="old data_* dir, or combo dir with met/tf/filter/")
    ap.add_argument("--rid", "--region-id", type=int, default=None,
                    help="only this region_id (default: all)")
    ap.add_argument("--lid", "--loc-id", type=int, default=None,
                    help="only this loc_id (default: all)")
    ap.add_argument("--kernel", default=None,
                    help="only this kernel name, e.g. bt (default: all)")
    ap.add_argument("--filt", default="", help="path to filt.x")
    ap.add_argument("--width", type=int, default=100)
    ap.add_argument("--nsamp", type=int, default=1000000)
    ap.add_argument("--plow", type=float, default=0.001)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    random.seed(args.seed)

    root = Path(args.data_root).resolve()
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    nspg = load_nspg(root)
    ngauge_tab = load_ngauge_table(root)
    filt_x = find_filt_x(root, args.filt or None)

    nsite = 0
    nmatch = 0
    if is_new_tree(root):
        jobs = collect_new_tree_jobs(
            root, ngauge_tab, nspg, args.rid, args.lid, args.kernel
        )
        nmatch = len(jobs)
        if not jobs and args.rid is None and args.lid is None:
            raise SystemExit(f"no tf/kernel.class_rx_lx/{{tfs,tfe}} under {root}")
        for k, c, met_dir, tfs_dir, tfe_dir, out_dir, rid, loc, ng in jobs:
            if run_one_site(
                k, c, met_dir, tfs_dir, tfe_dir, out_dir, rid, loc, ng,
                nspg, filt_x, args.width, args.nsamp, args.plow,
            ):
                nsite += 1
    else:
        met_dirs = find_kernel_dirs(root)
        tf_dirs = find_tf_dirs(root)
        pairs = []
        for (k, c), met_dir in met_dirs.items():
            tfs_dir = tf_dirs.get((k, c, "tfs"))
            tfe_dir = tf_dirs.get((k, c, "tfe"))
            if tfs_dir and tfe_dir:
                pairs.append((k, c, met_dir, tfs_dir, tfe_dir))
        if not pairs:
            raise SystemExit(f"no Kernel.CLASS + _tfs + _tfe trio under {root}")
        for k, c, met_dir, tfs_dir, tfe_dir in pairs:
            sites = collect_sites(
                k, c, tfs_dir, ngauge_tab, nspg, args.rid, args.lid
            )
            nmatch += len(sites)
            for rid, loc, ng in sites:
                out_dir = root / f"{met_dir.name}_filt" / f"r{rid}_l{loc}"
                if run_one_site(
                    k, c, met_dir, tfs_dir, tfe_dir, out_dir, rid, loc, ng,
                    nspg, filt_x, args.width, args.nsamp, args.plow,
                ):
                    nsite += 1
    if args.rid is not None or args.lid is not None:
        sel = []
        if args.rid is not None:
            sel.append(f"rid={args.rid}")
        if args.lid is not None:
            sel.append(f"lid={args.lid}")
        if nmatch == 0:
            raise SystemExit(
                f"no sites match {', '.join(sel)} under {root}"
            )
    if nsite == 0:
        raise SystemExit("no sites filtered")
    print(f"filtered {nsite} site(s) under {root}")


if __name__ == "__main__":
    main()
