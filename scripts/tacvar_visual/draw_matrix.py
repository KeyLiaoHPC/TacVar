"""draw_matrix.py — 为 run_matrix.sh 的 manifest 批量产出打印质量图（本机运行）。

用法:
  python3 draw_matrix.py <suite目录> [--kernel cg] [--skip-figures] [--only-timer-compare]

读取 <suite>/matrix/manifest.csv，对每个完成的组合:
  1) tr_vs_met / tf_cdf / gauge_overhead 单站点图 → <data_dir>/images/
  2) 每个 kernel 每个 (region_id, loc_id) 的 timer 对比图（同站点的
     tr CDF 叠加: 6 种 timer×counter 组合）→ <suite>/matrix/figures/

依赖: tacvar_visual.draw_print（数据 NFS 共享，形态同数据目录）。
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

from tacvar_visual.draw_print import (
    ACADEMIC_RC,
    REGION_NAMES,
    FIG_W_IN,
    FIG_H_IN,
    DPI,
    _auto_unit,
    load_tr_hist,
)

_COMBO_LABEL = {
    "T1": "clock_gettime",
    "T2": "clock_gettime + 4xPAPI",
    "T3": "cntvct_el0",
    "T4": "cntvct_el0 + 4xPAPI",
    "T5": "papi_get_real_nsec",
    "T6": "papi_get_real_nsec + 4xPAPI",
}
_COMBO_COLOR = {
    "T1": "#1f77b4",
    "T2": "#1f77b4",
    "T3": "#d62728",
    "T4": "#d62728",
    "T5": "#2ca02c",
    "T6": "#2ca02c",
}
_COMBO_LS = {"T1": "-", "T2": "--", "T3": "-", "T4": "--", "T5": "-", "T6": "--"}


def _tr_cdf(site_dir: Path) -> tuple[np.ndarray, np.ndarray] | None:
    tr_path = site_dir / "tr_hist.csv"
    if not tr_path.is_file():
        return None
    t, p = load_tr_hist(tr_path)
    if t.size < 2:
        return None
    edges = np.concatenate([t, [t[-1] + (t[-1] - t[0])]])
    cdf = np.cumsum(p)
    return edges, cdf


def _median_from_cdf(edges: np.ndarray, cdf: np.ndarray) -> float:
    return float(edges[np.searchsorted(cdf, 0.5)])


def load_manifest(suite: Path) -> list[dict]:
    mf = suite / "matrix" / "manifest.csv"
    if not mf.is_file():
        raise SystemExit(f"missing manifest: {mf}")
    rows = []
    with mf.open(newline="") as fp:
        for rec in csv.DictReader(fp):
            if rec.get("verification", "").startswith("FAIL"):
                continue
            rows.append(rec)
    return rows


def draw_timer_compare(
    suite: Path,
    rows: list[dict],
    *,
    kernel_filter: str = "",
    formats: tuple[str, ...] = ("pdf", "png"),
    dpi: int = DPI,
) -> int:
    """每个 kernel 每个 (rid,loc): 将已完成组合的 tr CDF 叠在一张 16:9 图。"""
    by_kernel: dict[str, list[dict]] = {}
    for r in rows:
        k = r["bench"].strip().upper()
        if kernel_filter and k != kernel_filter.upper():
            continue
        by_kernel.setdefault(k, []).append(r)

    out_root = suite / "matrix" / "figures"
    out_root.mkdir(parents=True, exist_ok=True)
    n_figs = 0

    for kernel, combos in sorted(by_kernel.items()):
        site_dict: dict[tuple[int, int], list[dict]] = {}
        for r in combos:
            data_dir = suite / r["data_dir"].strip()
            if not data_dir.is_dir():
                continue
            # FilT 输出目录名沿用数据目录的小写 Kernel.Class（如 bt.C_filt），
            # 不能按大写资源名拼接；扫描实际的 *_filt 目录。
            filt_roots = [
                d for d in data_dir.iterdir()
                if d.is_dir() and re.fullmatch(r"[A-Za-z0-9]+\.[A-Za-z]_filt", d.name)
            ]
            for filt_root in filt_roots:
                for sd in filt_root.iterdir():
                    if not sd.is_dir() or not sd.name.startswith("r") or "_l" not in sd.name:
                        continue
                    try:
                        rid, loc = sd.name[1:].split("_l")
                        key = (int(rid), int(loc))
                    except ValueError:
                        continue
                    site_dict.setdefault(key, []).append((r, sd))

        for (rid, loc), entries in sorted(site_dict.items()):
            if len({e[0]["tag"] for e in entries}) < 2:
                continue  # 至少两个 timer 组合才成对比图
            fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), dpi=dpi)
            for spine in ax.spines.values():
                spine.set_color("black")
                spine.set_linewidth(0.9)

            all_vals: list[float] = []
            plotted = []
            for r, sd in sorted(entries, key=lambda e: e[0]["tag"]):
                res = _tr_cdf(sd)
                if res is None:
                    continue
                edges, cdf = res
                all_vals.extend(edges.tolist())
                med = _median_from_cdf(edges, cdf)
                ax.plot(
                    edges / 1.0,
                    np.concatenate([[0.0], cdf]),
                    color=_COMBO_COLOR.get(r["tag"], "#333333"),
                    linestyle=_COMBO_LS.get(r["tag"], "-"),
                    linewidth=1.6,
                    drawstyle="steps-post",
                    label=f"{_COMBO_LABEL.get(r['tag'], r['tag'])}  (med {med / 1e6:.3g} ms)",
                )
                plotted.append(r["tag"])
            if not plotted:
                plt.close(fig)
                continue

            scale, unit = _auto_unit(np.asarray(all_vals, dtype=float))
            # 重新按单位缩放绘制（上面用 ns 画的，这里统一换算）
            ax.clear()
            for spine in ax.spines.values():
                spine.set_color("black")
                spine.set_linewidth(0.9)
            for r, sd in sorted(entries, key=lambda e: e[0]["tag"]):
                res = _tr_cdf(sd)
                if res is None:
                    continue
                edges, cdf = res
                med = _median_from_cdf(edges, cdf)
                ax.plot(
                    edges / scale,
                    np.concatenate([[0.0], cdf]),
                    color=_COMBO_COLOR.get(r["tag"], "#333333"),
                    linestyle=_COMBO_LS.get(r["tag"], "-"),
                    linewidth=1.6,
                    drawstyle="steps-post",
                    label=f"{_COMBO_LABEL.get(r['tag'], r['tag'])}  (med {med / scale:.3g} {unit})",
                )
            ax.axhline(0.5, color="0.5", linestyle=":", linewidth=0.9)
            rname = REGION_NAMES.get(rid, f"region {rid}")
            ax.set_xlabel(f"filtered real time tr ({unit})")
            ax.set_ylabel("cumulative probability")
            ax.set_ylim(-0.02, 1.05)
            ax.set_title(
                f"{kernel}.{entries[0][0]['class'].strip().upper()} · {rname} "
                f"(region {rid}, loc {loc}) — tr CDF per timer/counter",
                fontsize=16,
            )
            ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=11)
            fig.tight_layout()
            base = out_root / f"{kernel}_r{rid}_l{loc}_timer_compare"
            for fmt in formats:
                fig.savefig(f"{base}.{fmt}", dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            n_figs += 1
            print(f"saved {base}.{'/'.join(formats)}  (timer combos: {','.join(plotted)})")
    return n_figs


def main() -> None:
    ap = argparse.ArgumentParser(description="batch draw for run_matrix manifest")
    ap.add_argument("suite", help="NPB3.4-MPI 套件目录（含 matrix/manifest.csv）")
    ap.add_argument("--kernel", default="", help="只画某 kernel（如 cg）")
    ap.add_argument("--formats", default="pdf,png")
    ap.add_argument("--no-sites", action="store_true", help="跳过单站点图，只出 timer 对比图")
    args = ap.parse_args()

    suite = Path(args.suite).resolve()
    sys.path.insert(0, str(suite.parents[2] / "scripts"))  # <repo>/scripts
    from tacvar_visual.draw_print import (  # noqa: F401  (re-import for module path)
        draw_print_gauge_overhead,
        draw_print_tf_cdf,
        draw_print_tr_vs_met,
    )

    rows = load_manifest(suite)
    if not rows:
        print("no completed rows in manifest")
        return
    print(f"manifest: {len(rows)} completed rows")

    fmts = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    with plt.rc_context(ACADEMIC_RC):
        if not args.no_sites:
            for r in rows:
                if args.kernel and r["bench"].upper() != args.kernel.upper():
                    continue
                data_dir = suite / r["data_dir"].strip()
                if not data_dir.is_dir():
                    print(f"skip (no dir) {r['data_dir']}")
                    continue
                k, c = r["bench"].strip(), r["class"].strip().upper()
                try:
                    draw_print_tr_vs_met(data_dir, k, c, formats=fmts, ref_note=f"{r['timer']}")
                    draw_print_tf_cdf(data_dir, k, c, formats=fmts)
                    draw_print_gauge_overhead(data_dir, k, c, formats=fmts)
                except Exception as e:  # noqa: BLE001
                    print(f"site figs failed for {r['data_dir']}: {e}", file=sys.stderr)
        n = draw_timer_compare(suite, rows, kernel_filter=args.kernel, formats=fmts)
        print(f"timer-compare figures: {n}")


if __name__ == "__main__":
    main()