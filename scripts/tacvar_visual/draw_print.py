"""Print-quality figures for NPB-MPI TacVar filtered results.

Reads one data_<stamp> directory produced by the in-situ sampling flow:

    data_<stamp>/
      nspg.txt                  ns per subtraction-gauge step
      median.csv                kernel,class,region_id,loc_id,median,ngauge
      cg.C/                     per-rank raw measurement CSVs
      cg.C_filt/r<rid>_l<loc>/  FilT outputs: met.csv, tf.csv, tf_cdf.csv,
                                tr_hist.csv, sim_cdf.csv, er.out, ep.out

Figures: 16:9 print quality (PDF + PNG), Cambria serif, full axis labels,
median / percentile annotations, no grid.

Author: Hermes / TacVar project (draw_print_* series).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# ---------------------------------------------------------------------------
# Style: academic 16:9 slides, Cambria first
# ---------------------------------------------------------------------------

ACADEMIC_RC = {
    "font.family": "serif",
    "font.serif": [
        "Cambria",
        "Liberation Serif",
        "DejaVu Serif",
        "Times New Roman",
        "serif",
    ],
    "font.size": 13,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "axes.linewidth": 0.9,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
}

FIG_W_IN, FIG_H_IN = 13.333, 7.5  # 16:9
DPI = 300

# Region catalogs (names_* in common/tacvar_npb.c)
REGION_NAMES = {
    "is": {0: "total", 1: "rcomp", 2: "rcomm", 3: "verify"},
    "cg": {1: "total", 2: "conjg", 3: "rcomm", 4: "ncomm"},
    "mg": {
        1: "total", 2: "init", 3: "psinv", 4: "resid",
        5: "rprj3", 6: "interp", 7: "norm2u3", 8: "comm3", 9: "rcomm",
    },
    "ep": {1: "total", 2: "gpairs", 3: "randn", 4: "rcomm"},
    "bt": {
        1: "total", 2: "i/o", 3: "rhs", 4: "xsolve", 5: "ysolve",
        6: "zsolve", 7: "bpack", 8: "exch", 9: "xcomm", 10: "ycomm",
        11: "zcomm", 12: "enorm", 13: "iov",
    },
    "sp": {
        1: "total", 2: "rhs", 3: "xsolve", 4: "ysolve", 5: "zsolve",
        6: "bpack", 7: "exch", 8: "xcomm", 9: "ycomm", 10: "zcomm",
    },
    "lu": {
        1: "total", 2: "rhs", 3: "blts", 4: "buts", 5: "jacld",
        6: "jacu", 7: "exch", 8: "lcomm", 9: "ucomm", 10: "rcomm",
    },
    "ft": {
        1: "total", 2: "setup", 3: "fft", 4: "evolve", 5: "checksum",
        6: "fftlow", 7: "fftcopy", 8: "transpose",
        9: "transxzloc", 10: "transxzglo", 11: "transxzfin",
        12: "transxyloc", 13: "transxyglo", 14: "transxyfin",
        15: "synch", 16: "init",
    },
}

_MET_BLUE = "#1f77b4"
_MET_GREY = "0.35"
_TR_COLOR = "#c00000"
_SIM_COLOR = "#2ca02c"


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def _read_col(path: Path) -> np.ndarray:
    """One integer per line (met.csv / tf.csv)."""
    vals = []
    with path.open() as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            vals.append(int(float(line.split(",")[0])))
    return np.asarray(vals, dtype=np.float64)


def load_tr_hist(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """tr_hist.csv: (t, p) rows; last row p=0 marks the right boundary."""
    ts, ps = [], []
    with path.open() as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 2:
                continue
            ts.append(float(parts[0]))
            ps.append(float(parts[1]))
    t = np.asarray(ts, dtype=np.float64)
    p = np.asarray(ps, dtype=np.float64)
    # Drop the right-boundary marker row if present.
    if p.size > 0 and p[-1] == 0.0 and t.size > 1:
        t, p = t[:-1], p[:-1]
    return t, p


def load_tf_cdf(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """tf_cdf.csv: percentile, tf_ns (1000 rows)."""
    qs, vs = [], []
    with path.open() as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("percentile"):
                continue
            parts = line.split(",")
            qs.append(float(parts[0]))
            vs.append(float(parts[1]))
    return np.asarray(qs), np.asarray(vs)


def load_scalar(path: Path) -> float:
    with path.open() as fp:
        return float(fp.read().strip())


def load_median_table(root: Path) -> dict[tuple[int, int], tuple[int, int]]:
    """(rid, loc) -> (median_ns, ngauge) from median.csv."""
    out: dict[tuple[int, int], tuple[int, int]] = {}
    path = root / "median.csv"
    if not path.is_file():
        return out
    with path.open(newline="") as fp:
        for rec in csv.DictReader(fp):
            try:
                rid = int(rec["region_id"])
                loc = int(rec["loc_id"])
                med = int(float(rec["median"]))
                ng = int(float(rec["ngauge"]))
            except (KeyError, ValueError):
                continue
            out[(rid, loc)] = (med, ng)
    return out


def load_nspg(root: Path) -> float:
    path = root / "nspg.txt"
    with path.open() as fp:
        return float(fp.read().split()[0])


def list_filt_sites(root: Path) -> list[tuple[int, int, Path]]:
    """Sorted (rid, loc, site_dir) from <root>/<K>.<C>_filt/r*_l*/."""
    sites = []
    for kd in root.iterdir():
        m = re.match(r"^([A-Za-z0-9]+)\.([A-Za-z])_filt$", kd.name)
        if not m:
            continue
        kernel, klass = m.group(1).upper(), m.group(2).upper()
        for sd in kd.iterdir():
            m2 = re.fullmatch(r"r(\d+)_l(\d+)", sd.name)
            if not m2 or not sd.is_dir():
                continue
            sites.append((kernel, klass, int(m2.group(1)), int(m2.group(2)), sd))
    return sorted(sites)


def _auto_unit(values_ns: np.ndarray) -> tuple[float, str]:
    """Pick (scale, unit) so ticks are comfortable: ns / us / ms / s."""
    vmax = float(np.nanmax(values_ns)) if values_ns.size else 1.0
    if vmax >= 1e9:
        return 1e9, "s"
    if vmax >= 1e6:
        return 1e6, "ms"
    if vmax >= 1e3:
        return 1e3, r"$\mu$s"
    return 1.0, "ns"


def _format_annot(median: float, p95: float, unit: str, scale: float) -> str:
    return f"median {median / scale:.3g} {unit}  ·  P95 {p95 / scale:.3g} {unit}"


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------

def _new_figure() -> tuple[Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), dpi=DPI)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.9)
    return fig, ax


def resolve_filt_dir(root: Path, kernel: str, klass: str) -> Path:
    """Find <root>/<kernel>.<class>_filt (case-insensitive kernel)."""
    ku, cu = kernel.upper(), klass.upper()
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        m = re.fullmatch(r"([A-Za-z0-9]+)\.([A-Za-z])_filt", d.name)
        if m and m.group(1).upper() == ku and m.group(2).upper() == cu:
            return d
    raise ValueError(f"no FilT directory for {ku}.{cu} under {root}")


def _iter_sites(
    filt_root: Path, rid: int | None = None, lid: int | None = None
) -> list[tuple[int, int, Path]]:
    sites = []
    for sd in filt_root.iterdir():
        if not sd.is_dir():
            continue
        m = re.fullmatch(r"r(\d+)_l(\d+)", sd.name)
        if not m:
            continue
        r, l = int(m.group(1)), int(m.group(2))
        if rid is not None and r != rid:
            continue
        if lid is not None and l != lid:
            continue
        sites.append((r, l, sd))
    return sorted(sites)


def _site_title(kernel: str, klass: str, rid: int, loc: int, tag: str) -> str:
    rname = REGION_NAMES.get(kernel.lower(), {}).get(rid, f"region {rid}")
    return (
        f"{kernel.upper()}.{klass.upper()} · {rname} "
        f"(region {rid}, loc {loc}) — {tag}"
    )


def draw_print_tr_vs_met(
    data_root: str | Path,
    kernel: str = "cg",
    klass: str = "C",
    out_dir: str | Path | None = None,
    *,
    rid: int | None = None,
    lid: int | None = None,
    formats: tuple[str, ...] = ("pdf", "png"),
    ref_note: str = "",
) -> list[Figure]:
    """Per-site panel: raw met vs filtered tr, histogram (left) and CDF (right).

    Saves data_<stamp>/images/{K}.{C}_tr_vs_met_r{rid}_l{loc}.{fmt}.
    """
    root = Path(data_root)
    filt_root = resolve_filt_dir(root, kernel, klass)
    if out_dir is None:
        out_dir = root / "images"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    med_table = load_median_table(root)
    figures: list[Figure] = []
    kernel_u, klass_u = kernel.upper(), klass.upper()

    sites = _iter_sites(filt_root, rid, lid)
    if not sites:
        raise ValueError(f"no FilT sites under {filt_root} (rid={rid} lid={lid})")

    for rid, loc, site in sites:
        met_path = site / "met.csv"
        tr_path = site / "tr_hist.csv"
        sim_path = site / "sim_cdf.csv"
        if not (met_path.is_file() and tr_path.is_file()):
            print(f"skip {site.name}: missing met.csv/tr_hist.csv", file=sys.stderr)
            continue

        met = _read_col(met_path)
        t_tr, p_tr = load_tr_hist(tr_path)
        if t_tr.size < 2 or p_tr.size < 2:
            print(f"skip {site.name}: tr_hist degenerate", file=sys.stderr)
            continue

        all_ns = np.concatenate([met, t_tr])
        scale, unit = _auto_unit(all_ns)
        met_med, met_p95 = float(np.median(met)), float(np.percentile(met, 95))
        # tr CDF from histogram (cumulative), including right boundary
        tr_edges = np.concatenate([t_tr, [t_tr[-1] + (t_tr[-1] - t_tr[0] if t_tr.size > 1 else 1.0)]])
        tr_cdf = np.cumsum(p_tr)
        tr_med = tr_edges[np.searchsorted(tr_cdf, 0.5)]
        tr_p95 = tr_edges[np.searchsorted(tr_cdf, 0.95)]
        (med_med, med_ng) = med_table.get((rid, loc), (np.nan, np.nan))

        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(FIG_W_IN, FIG_H_IN), dpi=DPI
        )
        for ax in (ax1, ax2):
            for spine in ax.spines.values():
                spine.set_color("black")
                spine.set_linewidth(0.9)

        # ---- left: PDF
        n_bins = int(np.clip(round(np.sqrt(met.size)), 20, 120))
        ax1.hist(
            met / scale,
            bins=n_bins,
            density=True,
            color=_MET_GREY,
            alpha=0.45,
            edgecolor="black",
            linewidth=0.5,
            label=f"raw met (n={met.size:,})",
        )
        # tr histogram: each bin drawn with its own width, density = p / width
        tr_widths = np.diff(tr_edges) / scale
        ax1.bar(
            t_tr / scale,
            p_tr / tr_widths,
            width=tr_widths,
            align="edge",
            color=_TR_COLOR,
            alpha=0.55,
            edgecolor=_TR_COLOR,
            linewidth=0.4,
            label="filtered tr",
        )
        ax1.set_xlabel(f"time ({unit})")
        ax1.set_ylabel("probability density")
        ax1.legend(loc="upper right", frameon=True, framealpha=0.9)
        ax1.set_title("Histogram: before (met) vs after (tr)")

        # ---- right: CDF
        met_sorted = np.sort(met / scale)
        ax2.plot(
            met_sorted,
            np.arange(1, met_sorted.size + 1) / met_sorted.size,
            color=_MET_GREY,
            linewidth=1.6,
            label="raw met",
        )
        ax2.plot(
            tr_edges / scale,
            np.concatenate([[0.0], tr_cdf]),
            color=_TR_COLOR,
            linewidth=1.8,
            drawstyle="steps-post",
            label="filtered tr",
        )
        ax2.axhline(0.5, color="0.5", linestyle=":", linewidth=0.9)
        ax2.set_xlabel(f"time ({unit})")
        ax2.set_ylabel("cumulative probability")
        ax2.set_ylim(-0.02, 1.05)
        ax2.legend(loc="lower right", frameon=True, framealpha=0.9)
        ax2.set_title("CDF: before (met) vs after (tr)")

        fig.suptitle(
            _site_title(kernel, klass, rid, loc, "in-situ timing-fluctuation filter"),
            fontsize=17,
            y=0.98,
        )
        axis_note = (
            f"tuned tr: median {tr_med / scale:.3g} {unit}, P95 {tr_p95 / scale:.3g} {unit}  ·  "
            f"raw met: median {met_med / scale:.3g} {unit}, P95 {met_p95 / scale:.3g} {unit}"
        )
        if not np.isnan(med_ng):
            axis_note += f"  ·  ngauge={int(med_ng):d}"
        if ref_note:
            axis_note += f"  ·  {ref_note}"
        fig.text(0.5, 0.005, axis_note, ha="center", fontsize=11, color="0.15")
        fig.tight_layout(rect=(0, 0.025, 1, 0.96))

        base = out_dir / f"{kernel_u}.{klass_u}_tr_vs_met_r{rid}_l{loc}"
        for fmt in formats:
            fig.savefig(f"{base}.{fmt}", dpi=DPI, bbox_inches="tight")
        figures.append(fig)
        print(f"saved {base}.{'/'.join(formats)}  (tr med {tr_med/scale:.3g} {unit})")
        plt.close(fig)
    return figures


def draw_print_tf_cdf(
    data_root: str | Path,
    kernel: str = "cg",
    klass: str = "C",
    out_dir: str | Path | None = None,
    *,
    rid: int | None = None,
    lid: int | None = None,
    formats: tuple[str, ...] = ("pdf", "png"),
    ref_note: str = "",
) -> list[Figure]:
    """Per-site CDF of the subtracted timing fluctuation tf.

    Prefers tf_cdf.csv; falls back to the empirical CDF of tf.csv.
    Saves data_<stamp>/images/{K}.{C}_tf_cdf_r{rid}_l{loc}.{fmt}.
    """
    root = Path(data_root)
    filt_root = resolve_filt_dir(root, kernel, klass)
    if out_dir is None:
        out_dir = root / "images"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    figures: list[Figure] = []
    kernel_u, klass_u = kernel.upper(), klass.upper()
    for rid, loc, site in _iter_sites(filt_root, rid, lid):
        cdf_path = site / "tf_cdf.csv"
        tf_path = site / "tf.csv"
        if cdf_path.is_file():
            qs, vs = load_tf_cdf(cdf_path)
            qs = qs / 100.0
        elif tf_path.is_file():
            tf = _read_col(tf_path)
            if tf.size < 2:
                continue
            vs = np.sort(tf)
            qs = np.arange(1, vs.size + 1, dtype=np.float64) / vs.size
        else:
            continue
        scale, unit = _auto_unit(vs)

        fig, ax = _new_figure()
        ax.plot(vs / scale, qs, color=_MET_BLUE, linewidth=2.0)
        ax.axhline(0.5, color="0.5", linestyle=":", linewidth=0.9)
        ax.axvline(0.0, color="0.2", linestyle="--", linewidth=0.9)
        ax.set_xlabel(f"subtracted timing fluctuation tf ({unit})")
        ax.set_ylabel("cumulative probability")
        ax.set_ylim(-0.02, 1.05)
        ax.set_title(_site_title(kernel, klass, rid, loc, "timing fluctuation CDF"))
        note = (
            f"tf = elapsed_tf − ngauge·nspg (nspg={load_nspg(root):.4f} ns/step)"
        )
        if ref_note:
            note += f"  ·  {ref_note}"
        fig.text(0.5, 0.01, note, ha="center", fontsize=11, color="0.15")
        fig.tight_layout(rect=(0, 0.03, 1, 0.97))

        base = out_dir / f"{kernel_u}.{klass_u}_tf_cdf_r{rid}_l{loc}"
        for fmt in formats:
            fig.savefig(f"{base}.{fmt}", dpi=DPI, bbox_inches="tight")
        figures.append(fig)
        fig.clf()
        plt.close(fig)
    return figures


def draw_print_gauge_overhead(
    data_root: str | Path,
    kernel: str = "cg",
    klass: str = "C",
    out_dir: str | Path | None = None,
    *,
    rid: int | None = None,
    lid: int | None = None,
    formats: tuple[str, ...] = ("pdf", "png"),
) -> list[Figure]:
    """One bar chart: gauge overhead tmet = ngauge·nspg vs median(tr) per site.

    Saves data_<stamp>/images/{K}.{C}_gauge_overhead.{fmt}.
    """
    root = Path(data_root)
    filt_root = resolve_filt_dir(root, kernel, klass)
    if out_dir is None:
        out_dir = root / "images"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    nspg = load_nspg(root)
    med_table = load_median_table(root)
    kernel_u, klass_u = kernel.upper(), klass.upper()
    rows = []
    for rid, loc, site in _iter_sites(filt_root, rid, lid):
        tr_path = site / "tr_hist.csv"
        if not tr_path.is_file():
            continue
        t_tr, p_tr = load_tr_hist(tr_path)
        if t_tr.size < 2:
            continue
        tr_edges = np.concatenate([t_tr, [t_tr[-1] + (t_tr[-1] - t_tr[0])]])
        tr_cdf = np.cumsum(p_tr)
        tr_med = tr_edges[np.searchsorted(tr_cdf, 0.5)]
        med, ng = med_table.get((rid, loc), (np.nan, np.nan))
        if np.isnan(med) or ng < 0 or tr_med <= 0:
            continue
        tmet = ng * nspg
        rows.append((rid, loc, tr_med, tmet, tmet / tr_med * 100.0, ng))

    if not rows:
        raise ValueError("no usable sites for gauge-overhead chart")

    labels = [f"r{r}_l{l}" for r, l, *_ in rows]
    fracs = np.asarray([x[4] for x in rows])

    fig, ax = _new_figure()
    x = np.arange(len(rows))
    ax.bar(x, fracs, color=_TR_COLOR, alpha=0.7, edgecolor="black", linewidth=0.6)
    ax.axhline(1.0, color="0.4", linestyle="--", linewidth=0.9)
    for xi, (_, _, _, _, f, ng) in zip(x, rows):
        ax.text(xi, f + 0.04 * max(fracs.max(), 1.0), f"{f:.2f}%", ha="center", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("gauge overhead tmet / median(tr) (%)")
    ax.set_title(
        f"{kernel_u}.{klass_u} · in-situ gauge overhead per site "
        f"(nspg = {nspg:.4f} ns/step)"
    )
    ax.set_ylim(0, max(fracs.max() * 1.25, 2.0))
    fig.tight_layout()

    base = out_dir / f"{kernel_u}.{klass_u}_gauge_overhead"
    for fmt in formats:
        fig.savefig(f"{base}.{fmt}", dpi=DPI, bbox_inches="tight")
    print(f"saved {base}.{'/'.join(formats)}  ({len(rows)} sites)")
    plt.close(fig)
    return [fig]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="TacVar print-quality figures")
    ap.add_argument("data_root")
    ap.add_argument("--kernel", default="cg")
    ap.add_argument("--class", dest="klass", default="C")
    ap.add_argument("--rid", "--region-id", type=int, default=None)
    ap.add_argument("--lid", "--loc-id", type=int, default=None)
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--formats", default="pdf,png")
    ap.add_argument("--note", default="")
    ap.add_argument("--only", default="", help="comma list: tr_vs_met,tf_cdf,overhead")
    args = ap.parse_args()
    out = Path(args.out_dir) if args.out_dir else Path(args.data_root) / "images"
    fmts = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    only = set(f.strip() for f in args.only.split(",") if f.strip()) or None
    kw = dict(rid=args.rid, lid=args.lid, formats=fmts)

    with plt.rc_context(ACADEMIC_RC):
        if only is None or "tr_vs_met" in only:
            draw_print_tr_vs_met(
                args.data_root, args.kernel, args.klass, out,
                ref_note=args.note, **kw,
            )
        if only is None or "tf_cdf" in only:
            draw_print_tf_cdf(
                args.data_root, args.kernel, args.klass, out,
                ref_note=args.note, **kw,
            )
        if only is None or "overhead" in only:
            draw_print_gauge_overhead(
                args.data_root, args.kernel, args.klass, out, **kw,
            )
    print(f"figures written to {out.resolve()}")


if __name__ == "__main__":
    main()