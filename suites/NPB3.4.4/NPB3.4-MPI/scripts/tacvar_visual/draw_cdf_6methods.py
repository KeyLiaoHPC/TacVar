#!/usr/bin/env python3
"""
@file draw_cdf_6methods.py
@brief 三节点 × 6 计时组合的过滤前后 CDF 对比图（每站一张，12 条曲线）。

用法:
  python3 tacvar_visual/draw_cdf_6methods.py <matrix_glob> --out-dir <dir> \\
      [--formats png,pdf] [--dpi 300] [--stats <out.csv>]

matrix_glob: 引号包裹的 manifest.csv 通配（例如
  "/mnt/keylabmain/nfs/hpckey/03-Project/TacVar_NPC_*/suites/.../matrix/ext_*/T*/manifest.csv"）

每张图: 横轴 = 百分位 0..99（p99 截断，避免长尾拉扁主体）; 纵轴 = 耗时（ns/µs/ms 自动）;
6 种方法各一对曲线: 过滤前 met 原始分布（灰）与过滤后 tr 还原分布（强调色）。
er 标注在方法图例里。退化站（tr_hist 全 nan, 如 ep r3_l1）只画 met 并标注。

样式沿用 TacVar 学术风格（Cambria、白底、无网格、外刻度）；文本只用英文
（Cambria 无中文字形）。
"""
import argparse
import csv
import glob
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ACADEMIC_RC = {
    "font.family": "serif",
    "font.serif": ["Cambria", "Liberation Serif", "DejaVu Serif"],
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "axes.grid": False,
    "axes.edgecolor": "#404040",
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.top": False,
    "ytick.right": False,
    "legend.frameon": False,
    "legend.fontsize": 7.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "savefig.dpi": 300,
}
plt.rcParams.update(ACADEMIC_RC)

COMBO_FULL = {
    "T1": "native + none",
    "T2": "native + papi_read",
    "T3": "tick + none",
    "T4": "tick + papi_read",
    "T5": "papi + none",
    "T6": "papi + papi_read",
}
ORDER = ["T1", "T2", "T3", "T4", "T5", "T6"]
MET_COLOR = "#9a9a9a"
TR_COLORS = {
    "T1": "#1f77b4", "T2": "#ff7f0e", "T3": "#2ca02c",
    "T4": "#d62728", "T5": "#9467bd", "T6": "#8c564b",
}
TR_ALPHA = 0.85
MET_LW = 1.1
TR_LW = 1.6


def load_flat_csv(path: Path) -> np.ndarray | None:
    """met.csv / tf.csv 的单列 int64 数组。"""
    if not path.is_file():
        return None
    vals = []
    with open(path, newline="") as fp:
        for ln in fp:
            ln = ln.strip()
            if not ln:
                continue
            try:
                vals.append(int(float(ln.split()[0])))
            except (ValueError, IndexError):
                continue
    if not vals:
        return None
    return np.asarray(vals, dtype=np.int64)


def load_tr_hist(path: Path):
    """tr_hist.csv -> (bin_left_ns[], prob[])。格式: 'bin_left_edge_ns, prob'。
    最后一行的 p=0 是右边界标记，需丢弃；退化输出（全 nan / 全 0）返回 None。"""
    if not path.is_file():
        return None
    edges, probs = [], []
    with open(path) as fp:
        for ln in fp:
            ln = ln.strip().replace(",", " ")
            parts = ln.split()
            if len(parts) < 2:
                continue
            try:
                e = float(parts[0])
                p = float(parts[1])
            except ValueError:
                continue
            edges.append(e)
            probs.append(p)
    if len(edges) < 2:
        return None
    if np.isnan(probs).any() or (np.asarray(probs) <= 0).all():
        return None
    edges = np.asarray(edges[:-1])
    probs = np.asarray(probs[:-1])
    if probs.sum() <= 0:
        return None
    if edges[0] < 0.0 or not np.isfinite(edges).all():
        return None  # 负边界/非有限直方图 = 退化输出（如 er 巨大的组合）
    probs = probs / probs.sum()
    return edges, probs


def tr_percentiles(edges, probs, pct_grid):
    cdf = np.cumsum(probs)
    out = []
    for p in pct_grid:
        t = p / 100.0
        idx = np.searchsorted(cdf, t)
        idx = min(idx, len(edges) - 1)
        out.append(edges[idx])
    return np.asarray(out)


def met_percentiles(vals: np.ndarray, pct_grid) -> np.ndarray:
    return np.percentile(vals, pct_grid)


def yfmt(v):
    if v >= 1e6:
        return f"{v / 1e6:.2f} ms"
    if v >= 1e3:
        return f"{v / 1e3:.1f} us"
    return f"{v:.0f} ns"


def unit_scale(span_ns):
    if span_ns >= 1e6:
        return 1e6, "ms"
    if span_ns >= 1e3:
        return 1e3, "us"
    return 1.0, "ns"


def parse_manifest(mf: Path):
    """返回 [(node, bench, site, combo, er, data_dir), ...]"""
    node = "unknown"
    m = re.search(r"TacVar_NPC_([A-Za-z0-9]+)/suites", str(mf))
    combo = "?"
    m2 = re.search(r"/T([1-6])/manifest\.csv$", str(mf))
    if m:
        node = m.group(1)
    if m2:
        combo = "T" + m2.group(1)
    rows = []
    with open(mf, newline="") as fp:
        rd = csv.reader(fp)
        header = next(rd, None)
        if not header:
            return rows
        for rec in rd:
            if len(rec) < 12:
                continue
            bench, site, status = rec[0], rec[1], rec[2]
            if status not in ("SUCCESSFUL", "MEASURED") or site == "site":
                continue
            if site.upper().endswith(".C") and len(rec) > 3 and rec[3] == "OK":
                continue  # met 行（bench,BT.C,...）跳过
            er = rec[9] if len(rec) > 9 else ""
            data_dir = rec[10]
            rows.append((node, bench, site, combo, er, data_dir))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix_glob")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--formats", default="png,pdf")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--stats", default="")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pct_grid = np.arange(0, 100, dtype=float)  # 0..99，p99 截断
    combos = {t: {} for t in ORDER}
    sites = {}   # (node, bench, site) -> {combo: (er, data_dir)}

    for mf in sorted(glob.glob(args.matrix_glob)):
        for node, bench, site, combo, er, data_dir in parse_manifest(Path(mf)):
            if combo not in combos:
                continue
            key = (node, bench, site)
            sites.setdefault(key, {})[combo] = (er, data_dir)

    stats_rows = []
    for (node, bench, site), by_combo in sorted(sites.items()):
        fig, ax = plt.subplots(figsize=(9.0, 5.4))
        rid, loc = site[1:].split("_l")
        title = f"{node}  {bench}.C  site r{rid}_l{loc}   (percentile <= p99)"
        ax.set_title(title)
        ax.set_xlabel("percentile")
        ymin, ymax = np.inf, -np.inf
        plotted = 0
        for combo in ORDER:
            if combo not in by_combo:
                continue
            er, data_dir = by_combo[combo]
            filt = Path(data_dir) / f"{bench}.C_filt" / f"r{rid}_l{loc}"
            met = load_flat_csv(filt / "met.csv")
            hist = load_tr_hist(filt / "tr_hist.csv")
            if met is None:
                continue
            y_met = met_percentiles(met, pct_grid)
            ax.plot(pct_grid, y_met, color=MET_COLOR, lw=MET_LW, alpha=0.7,
                    label=f"{combo} met (er={er[:8] if er else '?'})" if hist is None
                    else f"{combo} met")
            plotted += 1
            if hist is not None:
                edges, probs = hist
                y_tr = tr_percentiles(edges, probs, pct_grid)
                ax.plot(pct_grid, y_tr, color=TR_COLORS[combo], lw=TR_LW,
                        alpha=TR_ALPHA, label=f"{combo} tr (er={er[:8] if er else '?'})")
                plotted += 1
            else:
                ax.annotate(f"{combo}: degenerate FilT output",
                            xy=(0.99, 0.02), xycoords="axes fraction",
                            ha="right", fontsize=7, color=TR_COLORS[combo])
            ymin = min(ymin, float(np.percentile(y_met, 1)))
            ymax = max(ymax, float(np.percentile(y_met, 99)))
            if hist is not None:
                ymax = max(ymax, float(np.percentile(y_tr, 99)))
        span = max(ymax - ymin, 1.0)
        scale, unit = unit_scale(span)
        ymin_s, ymax_s = ymin / scale, ymax / scale
        ax.set_ylim(max(0, ymin_s - 0.05 * span / scale), ymax_s + 0.05 * span / scale)
        ax.set_ylabel(f"time ({unit})")
        ax.legend(loc="upper left", ncol=2)
        fig.tight_layout()
        for fmt in args.formats.split(","):
            fig.savefig(out / f"{node}_{bench}_r{rid}_l{loc}.{fmt}", dpi=args.dpi)
        plt.close(fig)
        # 统计行
        for combo in ORDER:
            if combo in by_combo:
                er, data_dir = by_combo[combo]
                filt = Path(data_dir) / f"{bench}.C_filt" / f"r{rid}_l{loc}"
                met = load_flat_csv(filt / "met.csv")
                hist = load_tr_hist(filt / "tr_hist.csv")
                row = {"node": node, "bench": bench, "site": site, "combo": combo,
                       "er": er}
                if met is not None:
                    row.update({"tmet_p50": int(np.percentile(met, 50)),
                                "tmet_p99": int(np.percentile(met, 99))})
                if hist is not None:
                    tr = tr_percentiles(*hist, np.array([50.0, 99.0]))
                    row.update({"treal_p50": int(tr[0]), "treal_p99": int(tr[1])})
                stats_rows.append(row)
        print(f"drew {node} {bench} {site}: {plotted} curves")

    if args.stats and stats_rows:
        with open(args.stats, "w", newline="") as fp:
            w = csv.DictWriter(fp, fieldnames=list(stats_rows[0].keys()))
            w.writeheader()
            w.writerows(stats_rows)
        print(f"wrote {args.stats} ({len(stats_rows)} rows)")


if __name__ == "__main__":
    main()