#!/usr/bin/env python3
"""
@file draw_cdf_percentile.py
@brief Per-site percentile-CDF figures (tmet raw vs treal filtered) + stats CSV.

Usage:
  python3 tacvar_visual/draw_cdf_percentile.py <data_root> --out-dir <dir> [--csv <path>]

<data_root> = 矩阵 data 目录（NPC_experiments/<stamp>/matrix/data），其下每个子目录是
  一棵 bench 数据树（含 <k>.<c>_filt/r<rid>_l<loc>/{met.csv,tr_hist.csv,...}）。

图：x 轴为百分位点 0..99（纵轴不被长尾拉高），y 轴为时间（ns/µs/ms 自动单位）；
   tmet（原始 met 采样）与 treal（FilT 估计真实分布）两条曲线同图对比。
统计：min/p50/p95（ns）写入 CSV。
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SITE_RE = re.compile(r"r(\d+)_l(\d+)$")
FONT_CANDIDATES = ["Cambria", "Liberation Serif", "DejaVu Serif", "Times New Roman"]
ACADEMIC_RC = {
    "font.family": "serif",
    "font.serif": FONT_CANDIDATES,
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "legend.fontsize": 13,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.grid": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.top": False,
    "ytick.right": False,
    "axes.linewidth": 1.0,
    "legend.frameon": False,
    "lines.linewidth": 1.8,
}

T_MET_COLOR = (0.45, 0.45, 0.45)      # 原始 met 灰色
T_REAL_COLOR = (0.72, 0.18, 0.12)     # FilT treal 强调色（深红）

UNIT_LABEL = [("ns", 1e0), ("µs", 1e3), ("ms", 1e6), ("s", 1e9)]


def load_tr_hist(path: Path):
    """tr_hist.csv: 'bin_left_edge_ns, prob'（每行时间带尾逗号）；末行 p=0 为右边界标记，
    丢弃。概率含 nan（FilT 退化输出，如 ep r3_l1）时返回 (None, None)。返回 (t_edges, prob)。"""
    t, p = [], []
    with path.open() as fp:
        for ln in fp:
            ln = ln.strip()
            if not ln:
                continue
            fields = ln.replace(",", " ").split()
            t.append(float(fields[0]))
            p.append(float(fields[1]))
    if t and p[-1] == 0.0:
        t = t[:-1]
        p = p[:-1]
    a = np.asarray(t, dtype=np.float64)
    b = np.asarray(p, dtype=np.float64)
    s = b.sum()
    if s > 0.0:
        b = b / s
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return None, None
    return a, b


def tr_quantiles(t, p, qs):
    """从直方图插值分位数（qs∈[0,1））。bin 左边界 t_i，视为连续直方图。"""
    c = np.cumsum(p)
    last_w = (t[-1] - t[-2]) if len(t) > 1 else max(t[-1] * 0.01, 1.0)
    edges = np.append(t, t[-1] + last_w)
    widths = np.diff(edges)
    out = []
    for q in qs:
        i = int(np.searchsorted(c, q))
        i = min(i, len(t) - 1)
        base = c[i] - p[i]
        frac = 0.0 if p[i] <= 0.0 else min(max((q - base) / p[i], 0.0), 1.0)
        out.append(t[i] + frac * widths[i])
    return np.asarray(out)


def auto_unit(v_ns: float) -> tuple[str, float]:
    for name, scale in UNIT_LABEL:
        if v_ns < scale * 1000.0 or name == "s":
            return name, scale
    return "s", 1e9


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", help="矩阵 data 目录（含各 bench 数据树）")
    ap.add_argument("--out-dir", default="figures/cdf_percentile")
    ap.add_argument("--csv", default="", help="统计 CSV 输出路径（默认 <data_root> 同级 tables/）")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--formats", default="png,pdf", help="逗号分隔: png,pdf")
    args = ap.parse_args()

    root = Path(args.data_root).resolve()
    if not root.is_dir():
        raise SystemExit(f"not a dir: {root}")
    out_dir = Path(args.out_dir).resolve() if args.out_dir.startswith("/") else root.parent / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = Path(args.csv) if args.csv else root.parent / "tables" / "tmet_treal_stats.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fmts = [f.strip() for f in args.formats.split(",") if f.strip()]

    plt.rcParams.update(ACADEMIC_RC)

    sites = []
    for bench_root in sorted(root.iterdir()):
        if not bench_root.is_dir():
            continue
        for filt in sorted(bench_root.glob("*_filt/r*_l*")):
            m = SITE_RE.search(filt.name)
            if not m:
                continue
            met_csv = filt / "met.csv"
            tr_csv = filt / "tr_hist.csv"
            if not (met_csv.is_file() and tr_csv.is_file()):
                print(f"skip incomplete {filt}", file=sys.stderr)
                continue
            met = np.loadtxt(met_csv, dtype=np.int64)
            t, p = load_tr_hist(tr_csv)
            if met.size == 0:
                print(f"skip empty {filt}", file=sys.stderr)
                continue
            tr_valid = t is not None and t.size > 0
            ids = bench_root.name.split("_")[0].split(".")[0]
            bench, klass = ids, (bench_root.name.split(".")[1] if "." in bench_root.name else "C")
            rid, loc = int(m.group(1)), int(m.group(2))

            qs = np.arange(0, 100)                 # 0..99 百分位
            q_met = np.percentile(met, qs)
            q_tr = tr_quantiles(t, p, qs / 100.0) if tr_valid else None
            if q_tr is not None and not np.isfinite(q_tr).all():
                # FilT 退化输出（nan 概率）：treal 不可统计，图只画 tmet
                tr_valid = False
                q_tr = None
            stats = dict(
                n_met=met.size, tmet_min=int(met.min()), tmet_p50=int(np.median(met)),
                tmet_p95=int(np.percentile(met, 95)),
                treal_min=float(t[0]) if tr_valid else float("nan"),
                treal_p50=float(tr_quantiles(t, p, [0.5])[0]) if tr_valid else float("nan"),
                treal_p95=float(tr_quantiles(t, p, [0.95])[0]) if tr_valid else float("nan"),
            )
            sites.append((bench, klass, rid, loc, dict(stats)))

            # ---- 图 ----
            ymax = max(q_met[-1], q_tr[-1] if q_tr is not None else q_met[-1])
            unit, scale = auto_unit(ymax)
            fig, ax = plt.subplots(figsize=(13.333, 7.5))
            ax.plot(qs, q_met / scale, color=T_MET_COLOR, lw=1.8, alpha=0.85,
                    label=f"tmet (raw, n={met.size:,})")
            if tr_valid:
                ax.plot(qs, q_tr / scale, color=T_REAL_COLOR, lw=2.2,
                        label="treal (FilT)")
            else:
                ax.text(0.98, 0.05, "treal invalid: degenerate FilT output (nan)",
                        transform=ax.transAxes, ha="right", va="bottom",
                        color=T_REAL_COLOR, fontsize=13)
            ax.set_xlim(0, 99)
            ax.set_xlabel("Percentile (%)")
            ax.set_ylabel(f"Time ({unit})")
            ax.set_title(f"{bench}.{klass} r{rid}_l{loc} — met vs filtered-tr (percentile CDF)")
            ax.legend(loc="upper left")
            ax.tick_params(axis="both", labelsize=13)
            base = out_dir / f"{bench}_r{rid}_l{loc}"
            for fmt in fmts:
                fig.savefig(f"{base}.{fmt}", dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)

    if not sites:
        raise SystemExit(f"no sites found under {root}")

    # ---- 统计 CSV + 终端表 ----
    with csv_path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["bench", "class", "region_id", "loc_id", "n_met",
                    "tmet_min_ns", "tmet_p50_ns", "tmet_p95_ns",
                    "treal_min_ns", "treal_p50_ns", "treal_p95_ns"])
        for bench, klass, rid, loc, s in sorted(sites):
            w.writerow([bench, klass, rid, loc, s["n_met"],
                        s["tmet_min"], s["tmet_p50"], s["tmet_p95"],
                        f"{s['treal_min']:.0f}", f"{s['treal_p50']:.0f}", f"{s['treal_p95']:.0f}"])
    print(f"wrote {len(sites)} rows -> {csv_path}")
    print(f"figures -> {out_dir} ({len(sites)} sites, formats={','.join(fmts)})")

    def fmt_ns(v):
        if v >= 1e6:
            return f"{v / 1e6:.3f} ms"
        if v >= 1e3:
            return f"{v / 1e3:.2f} µs"
        return f"{v:.0f} ns"

    print(f"\n{'site':16s} {'n_met':>9s} {'tmet min/p50/p95':>24s} {'treal min/p50/p95':>24s}")
    for bench, klass, rid, loc, s in sorted(sites):
        tm = f"{fmt_ns(s['tmet_min'])} / {fmt_ns(s['tmet_p50'])} / {fmt_ns(s['tmet_p95'])}"
        if np.isfinite(s["treal_p50"]):
            tr = f"{fmt_ns(s['treal_min'])} / {fmt_ns(s['treal_p50'])} / {fmt_ns(s['treal_p95'])}"
        else:
            tr = "无效 (FilT 退化)"
        print(f"{bench + ' r' + str(rid) + '_l' + str(loc):16s} {s['n_met']:>9,d} {tm:>24s} {tr:>24s}")


if __name__ == "__main__":
    main()