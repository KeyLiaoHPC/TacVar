"""Pooled histogram, PDF, and CDF plots for NPB-MPI TacVar measurements."""

from __future__ import annotations

from typing import Iterable, List, NamedTuple, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from .heatmap import (
    _ACADEMIC_RC,
    _DEFAULT_BARS,
    _build_rank_matrix,
    _display_notebook_figures,
    _format_edge_label,
    _format_host_label,
    _green_white_red_cmap,
    _ratio_edges_from_bars,
    _ratio_from_refs,
    _resolve_xrange,
)
from .io import (
    load_measurement_frame,
    load_timer_info_table,
    resolve_kernel_dir,
    select_ids,
    short_hosts_from_frame,
)

_FIG_H_IN = 10.24
_FIG_W_IN = 10.24
_FIG_DPI = 100


class _DistView(NamedTuple):
    rid: int
    loc: int
    rname: str
    host_label: str
    kernel_name: str
    kernel_class: str
    fom: str
    n_ranks: int
    n_pts: int
    n_full_cols: int
    n_samples: int
    ratio_flat: np.ndarray


def _style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)


def _new_figure():
    fig, ax = plt.subplots(figsize=(_FIG_W_IN, _FIG_H_IN), dpi=_FIG_DPI)
    fig.set_dpi(_FIG_DPI)
    _style_axes(ax)
    return fig, ax


def _title(view: _DistView, kind: str) -> str:
    return (
        f"{view.fom} {kind} of {view.kernel_name}.{view.kernel_class} "
        f"on {view.host_label} (region={view.rname}, loc_id={view.loc})"
    )


def _pct_from_ratio(ratio_flat: np.ndarray) -> np.ndarray:
    return (ratio_flat - 1.0) * 100.0


def _pct_xlim(pct: np.ndarray, bars: Sequence[float]) -> tuple:
    lo = float(np.min(pct))
    hi = float(np.max(pct))
    span = hi - lo
    pad = 1.0 if span <= 0.0 else max(1.0, 0.05 * span)
    bars_arr = np.asarray(bars, dtype=float)
    x0 = min(lo - pad, float(bars_arr[0]))
    x1 = max(hi + pad, float(bars_arr[-1]))
    if x1 <= x0:
        x1 = x0 + 2.0
    return (x0, x1)


def _paint_gwr_guides(ax, bars: Sequence[float], cmap, n_bins: int, xlim) -> None:
    edges_pct = np.concatenate(
        [[-np.inf], np.asarray(bars, dtype=float), [np.inf]]
    )
    x0, x1 = xlim
    for i in range(n_bins):
        lo = edges_pct[i]
        hi = edges_pct[i + 1]
        if np.isneginf(lo):
            lo = x0
        if np.isposinf(hi):
            hi = x1
        if hi <= x0 or lo >= x1:
            continue
        ax.axvspan(
            max(float(lo), x0),
            min(float(hi), x1),
            color=cmap(i),
            alpha=0.18,
            lw=0,
            zorder=0,
        )
    for p in bars:
        ax.axvline(
            float(p), color="0.4", linestyle="--", linewidth=0.6, zorder=2
        )
    ax.axvline(0.0, color="0.2", linestyle="-", linewidth=0.7, zorder=2)


def _collect_dist_views(
    data_root: str,
    kernel_name: str,
    kernel_class: str,
    *,
    hosts: Optional[Sequence[str]],
    region_ids: Optional[Iterable[int]],
    loc_ids: Optional[Iterable[int]],
    ranks: Optional[Iterable[int]],
    fom: str,
    xrange: Union[int, list, tuple],
    ref_section_base: bool,
    ref_key: str,
    ref_base: str,
    bars: Optional[Sequence[float]],
) -> tuple:
    if hosts is None:
        hosts = ["all"]
    if region_ids is None:
        region_ids = range(0)
    if loc_ids is None:
        loc_ids = range(0)
    if ranks is None:
        ranks = range(0)
    if bars is None:
        bars = list(_DEFAULT_BARS)

    kernel_dir = resolve_kernel_dir(data_root, kernel_name, kernel_class)
    kernel_label = f"{kernel_name}.{kernel_class}"
    timer_table = load_timer_info_table(kernel_dir)
    print(f"timer_info for {kernel_label} ({kernel_dir})")
    if timer_table.empty:
        print("  (timer_info.csv missing or empty)")
        region_names = {}
    else:
        print(timer_table.to_string(index=False))
        region_names = {}
        if "region_id" in timer_table.columns and "name" in timer_table.columns:
            region_names = {
                int(r): str(n)
                for r, n in zip(timer_table["region_id"], timer_table["name"])
            }

    df = load_measurement_frame(kernel_dir, hosts=hosts)
    for col in ("region_id", "loc_id", "rank", fom):
        if col not in df.columns:
            raise ValueError(
                f"required column {col!r} missing from CSVs under {kernel_dir}"
            )

    selected_rids = select_ids(df, "region_id", region_ids)
    selected_locs = select_ids(df, "loc_id", loc_ids)
    selected_ranks = select_ids(df, "rank", ranks)
    host_label = _format_host_label(short_hosts_from_frame(df))
    rank_set = set(selected_ranks)

    views: List[_DistView] = []
    for rid in selected_rids:
        rname = region_names.get(rid, f"r{rid}")
        for loc in selected_locs:
            sub = df[(df["region_id"] == rid) & (df["loc_id"] == loc)]
            if sub.empty:
                continue
            sub = sub[sub["rank"].isin(rank_set)]
            if sub.empty:
                continue
            mat, rank_arr = _build_rank_matrix(sub, fom)
            col_idx = _resolve_xrange(int(mat.shape[1]), xrange)
            draw = mat[:, col_idx]
            n_ranks = len(rank_arr)
            n_pts = int(draw.shape[1])
            n_samples = int(np.isfinite(draw).sum())
            print(
                f"region_id={rid} ({rname}), loc_id={loc}: "
                f"{n_ranks} ranks, {n_pts} data points"
                + (
                    f" of {mat.shape[1]}"
                    if n_pts != int(mat.shape[1])
                    else ""
                )
                + (
                    f" ({n_samples} samples)"
                    if n_samples != n_ranks * n_pts
                    else ""
                )
            )
            ratio = _ratio_from_refs(
                mat, col_idx, ref_key, ref_base, bool(ref_section_base)
            )
            finite = ratio[np.isfinite(ratio)]
            if finite.size == 0:
                raise ValueError(
                    f"no finite value/ref samples for region_id={rid}, "
                    f"loc_id={loc} after filtering"
                )
            views.append(
                _DistView(
                    rid=rid,
                    loc=loc,
                    rname=rname,
                    host_label=host_label,
                    kernel_name=kernel_name,
                    kernel_class=kernel_class,
                    fom=fom,
                    n_ranks=n_ranks,
                    n_pts=n_pts,
                    n_full_cols=int(mat.shape[1]),
                    n_samples=int(finite.size),
                    ratio_flat=np.asarray(finite, dtype=float),
                )
            )

    if not views:
        raise ValueError(
            "no (region_id, loc_id) pairs to plot after filtering "
            f"(region_ids={selected_rids}, loc_ids={selected_locs}, "
            f"ranks={selected_ranks})"
        )
    return views, list(bars)


def draw_histogram_npb_mpi(
    data_root: str = "",
    kernel_name: str = "",
    kernel_class: str = "",
    *,
    hosts: Optional[Sequence[str]] = None,
    region_ids: Optional[Iterable[int]] = None,
    loc_ids: Optional[Iterable[int]] = None,
    ranks: Optional[Iterable[int]] = None,
    fom: str = "elapsed_ns",
    xrange: Union[int, list, tuple] = -1,
    ref_section_base: bool = True,
    ref_key: str = "median",
    ref_base: str = "all",
    bars: Optional[Sequence[float]] = None,
) -> List[Figure]:
    """
    Draw discrete green-white-red histograms of value/ref for NPB-MPI CSVs.

    One figure per (region_id, loc_id). Selected ranks are pooled into a
    single count histogram. X = percent-vs-ref bins from ``bars``; Y = count.
    Bar color is the same discrete green-white-red map as the heatmap.

    Parameters
    ----------
    data_root, kernel_name, kernel_class
        Path pieces for ``{data_root}/{kernel_name}.{kernel_class}/``.
        Remaining arguments are keyword-only.
    hosts
        ``['all']`` / empty / None for all hosts; else short_host prefixes.
    region_ids, loc_ids
        Ids to plot; empty means all present (intersection that exists in data).
    ranks
        MPI ranks to include; empty means all ranks present. Ratios from
        selected ranks are flattened into one sample set (no per-rank series).
    fom
        Numeric CSV column (default ``elapsed_ns``).
    xrange
        Columns to include. ``-1`` (default) is all samples. A ``list`` of ints
        selects those indices. A ``tuple (a, b)`` of ints is the half-open
        slice ``[a, b)``; of floats in ``[0, 1]`` it is
        ``[int(n*a), int(n*b))``. Integer ``(0, 1)`` is only column 0.
    ref_section_base
        If True, compute refs on the ``xrange`` window; if False, on the
        full rank-by-sample matrix. Display still uses ``xrange``.
    ref_key
        ``max``, ``min``, ``median``, ``average``, or ``pxx`` percentile.
    ref_base
        Sample set for ``ref_key``: ``all`` (one ref), ``x`` (all ranks of
        this seq), or ``y`` (seqs of this rank, scoped by
        ``ref_section_base``). Applied per cell, then pooled.
    bars
        Strictly increasing percent-deviation cut points. Number of bins is
        ``len(bars)+1`` after wrapping with ``+/- inf``. The bin that
        contains ratio 1.0 is white.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    views, bars_used = _collect_dist_views(
        data_root,
        kernel_name,
        kernel_class,
        hosts=hosts,
        region_ids=region_ids,
        loc_ids=loc_ids,
        ranks=ranks,
        fom=fom,
        xrange=xrange,
        ref_section_base=ref_section_base,
        ref_key=ref_key,
        ref_base=ref_base,
        bars=bars,
    )
    edges, white_idx = _ratio_edges_from_bars(bars_used)
    n_bins = len(edges) - 1
    cmap = _green_white_red_cmap(n_bins, white_idx)
    cut_points = edges[1:-1]
    colors = [cmap(i) for i in range(n_bins)]
    tick_pos = np.arange(n_bins + 1) - 0.5
    tick_labels = [_format_edge_label(e) for e in edges]

    figures: List[Figure] = []
    with plt.rc_context(_ACADEMIC_RC):
        for view in views:
            bin_idx = np.digitize(view.ratio_flat, cut_points, right=True)
            counts = np.bincount(bin_idx, minlength=n_bins)
            fig, ax = _new_figure()
            ax.bar(
                np.arange(n_bins),
                counts,
                width=0.9,
                color=colors,
                edgecolor="black",
                linewidth=0.6,
                zorder=3,
            )
            ax.set_xlim(-0.6, n_bins - 0.4)
            ax.set_xticks(tick_pos)
            ax.set_xticklabels(tick_labels)
            ax.set_xlabel("percent vs ref")
            ax.set_ylabel("count")
            ax.set_title(_title(view, "histogram"))
            fig.tight_layout()
            figures.append(fig)

    _display_notebook_figures(figures)
    return figures


def draw_pdf_npb_mpi(
    data_root: str = "",
    kernel_name: str = "",
    kernel_class: str = "",
    *,
    hosts: Optional[Sequence[str]] = None,
    region_ids: Optional[Iterable[int]] = None,
    loc_ids: Optional[Iterable[int]] = None,
    ranks: Optional[Iterable[int]] = None,
    fom: str = "elapsed_ns",
    xrange: Union[int, list, tuple] = -1,
    ref_section_base: bool = True,
    ref_key: str = "median",
    ref_base: str = "all",
    bars: Optional[Sequence[float]] = None,
) -> List[Figure]:
    """
    Draw a pooled PDF of percent-vs-ref for NPB-MPI TacVar CSVs.

    One figure per (region_id, loc_id). Selected ranks are pooled into a
    single density histogram of ``(value/ref - 1) * 100``. ``bars`` are
    vertical guides with light green-white-red bands (same bins as the
    heatmap). No per-rank series.

    Parameters
    ----------
    data_root, kernel_name, kernel_class
        Path pieces for ``{data_root}/{kernel_name}.{kernel_class}/``.
        Remaining arguments are keyword-only.
    hosts
        ``['all']`` / empty / None for all hosts; else short_host prefixes.
    region_ids, loc_ids
        Ids to plot; empty means all present (intersection that exists in data).
    ranks
        MPI ranks to include; empty means all ranks present. Ratios from
        selected ranks are flattened into one sample set (no per-rank series).
    fom
        Numeric CSV column (default ``elapsed_ns``).
    xrange
        Columns to include. ``-1`` (default) is all samples. A ``list`` of ints
        selects those indices. A ``tuple (a, b)`` of ints is the half-open
        slice ``[a, b)``; of floats in ``[0, 1]`` it is
        ``[int(n*a), int(n*b))``. Integer ``(0, 1)`` is only column 0.
    ref_section_base
        If True, compute refs on the ``xrange`` window; if False, on the
        full rank-by-sample matrix. Display still uses ``xrange``.
    ref_key
        ``max``, ``min``, ``median``, ``average``, or ``pxx`` percentile.
    ref_base
        Sample set for ``ref_key``: ``all`` (one ref), ``x`` (all ranks of
        this seq), or ``y`` (seqs of this rank, scoped by
        ``ref_section_base``). Applied per cell, then pooled.
    bars
        Strictly increasing percent-deviation cut points. Used as vertical
        guides and GWR background bands, not as PDF histogram edges.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    views, bars_used = _collect_dist_views(
        data_root,
        kernel_name,
        kernel_class,
        hosts=hosts,
        region_ids=region_ids,
        loc_ids=loc_ids,
        ranks=ranks,
        fom=fom,
        xrange=xrange,
        ref_section_base=ref_section_base,
        ref_key=ref_key,
        ref_base=ref_base,
        bars=bars,
    )
    edges, white_idx = _ratio_edges_from_bars(bars_used)
    n_bins = len(edges) - 1
    cmap = _green_white_red_cmap(n_bins, white_idx)

    figures: List[Figure] = []
    with plt.rc_context(_ACADEMIC_RC):
        for view in views:
            pct = _pct_from_ratio(view.ratio_flat)
            xlim = _pct_xlim(pct, bars_used)
            n_pdf_bins = int(np.clip(round(np.sqrt(pct.size)), 16, 80))
            dens, pdf_edges = np.histogram(
                pct, bins=n_pdf_bins, range=xlim, density=True
            )
            fig, ax = _new_figure()
            _paint_gwr_guides(ax, bars_used, cmap, n_bins, xlim)
            ax.stairs(
                dens,
                pdf_edges,
                fill=True,
                color="0.25",
                alpha=0.45,
                linewidth=0.0,
                zorder=3,
            )
            ax.stairs(
                dens,
                pdf_edges,
                fill=False,
                color="black",
                linewidth=1.1,
                zorder=4,
            )
            ax.set_xlim(xlim)
            ax.set_xlabel("percent vs ref")
            ax.set_ylabel("probability density")
            ax.set_title(_title(view, "pdf"))
            fig.tight_layout()
            figures.append(fig)

    _display_notebook_figures(figures)
    return figures


def draw_cdf_npb_mpi(
    data_root: str = "",
    kernel_name: str = "",
    kernel_class: str = "",
    *,
    hosts: Optional[Sequence[str]] = None,
    region_ids: Optional[Iterable[int]] = None,
    loc_ids: Optional[Iterable[int]] = None,
    ranks: Optional[Iterable[int]] = None,
    fom: str = "elapsed_ns",
    xrange: Union[int, list, tuple] = -1,
    ref_section_base: bool = True,
    ref_key: str = "median",
    ref_base: str = "all",
    bars: Optional[Sequence[float]] = None,
) -> List[Figure]:
    """
    Draw a pooled empirical CDF of percent-vs-ref for NPB-MPI TacVar CSVs.

    One figure per (region_id, loc_id). Selected ranks are pooled into a
    single ECDF of ``(value/ref - 1) * 100``. ``bars`` are vertical guides
    with light green-white-red bands (same bins as the heatmap). A dotted
    line marks cumulative probability 0.5. No per-rank series.

    Parameters
    ----------
    data_root, kernel_name, kernel_class
        Path pieces for ``{data_root}/{kernel_name}.{kernel_class}/``.
        Remaining arguments are keyword-only.
    hosts
        ``['all']`` / empty / None for all hosts; else short_host prefixes.
    region_ids, loc_ids
        Ids to plot; empty means all present (intersection that exists in data).
    ranks
        MPI ranks to include; empty means all ranks present. Ratios from
        selected ranks are flattened into one sample set (no per-rank series).
    fom
        Numeric CSV column (default ``elapsed_ns``).
    xrange
        Columns to include. ``-1`` (default) is all samples. A ``list`` of ints
        selects those indices. A ``tuple (a, b)`` of ints is the half-open
        slice ``[a, b)``; of floats in ``[0, 1]`` it is
        ``[int(n*a), int(n*b))``. Integer ``(0, 1)`` is only column 0.
    ref_section_base
        If True, compute refs on the ``xrange`` window; if False, on the
        full rank-by-sample matrix. Display still uses ``xrange``.
    ref_key
        ``max``, ``min``, ``median``, ``average``, or ``pxx`` percentile.
    ref_base
        Sample set for ``ref_key``: ``all`` (one ref), ``x`` (all ranks of
        this seq), or ``y`` (seqs of this rank, scoped by
        ``ref_section_base``). Applied per cell, then pooled.
    bars
        Strictly increasing percent-deviation cut points. Used as vertical
        guides and GWR background bands.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    views, bars_used = _collect_dist_views(
        data_root,
        kernel_name,
        kernel_class,
        hosts=hosts,
        region_ids=region_ids,
        loc_ids=loc_ids,
        ranks=ranks,
        fom=fom,
        xrange=xrange,
        ref_section_base=ref_section_base,
        ref_key=ref_key,
        ref_base=ref_base,
        bars=bars,
    )
    edges, white_idx = _ratio_edges_from_bars(bars_used)
    n_bins = len(edges) - 1
    cmap = _green_white_red_cmap(n_bins, white_idx)

    figures: List[Figure] = []
    with plt.rc_context(_ACADEMIC_RC):
        for view in views:
            pct = _pct_from_ratio(view.ratio_flat)
            xs = np.sort(pct)
            ys = np.arange(1, xs.size + 1, dtype=float) / float(xs.size)
            xlim = _pct_xlim(pct, bars_used)
            fig, ax = _new_figure()
            _paint_gwr_guides(ax, bars_used, cmap, n_bins, xlim)
            ax.step(xs, ys, where="post", color="black", linewidth=1.2, zorder=4)
            ax.axhline(0.5, color="0.5", linestyle=":", linewidth=0.8, zorder=3)
            ax.set_xlim(xlim)
            ax.set_ylim(0.0, 1.02)
            ax.set_xlabel("percent vs ref")
            ax.set_ylabel("cumulative probability")
            ax.set_title(_title(view, "cdf"))
            fig.tight_layout()
            figures.append(fig)

    _display_notebook_figures(figures)
    return figures
