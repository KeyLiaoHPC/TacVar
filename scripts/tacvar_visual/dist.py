"""Pooled histogram, PDF, and CDF plots for NPB-MPI TacVar measurements."""

from __future__ import annotations

import os
from collections import OrderedDict
from typing import Iterable, List, NamedTuple, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from .heatmap import (
    _ACADEMIC_RC,
    _DEFAULT_BARS,
    _build_rank_matrix,
    _compute_ref,
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
    normalize_data_roots,
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
    value_flat: np.ndarray
    ref: float
    run_label: str


def _style_axes(ax) -> None:
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)


def _new_figure():
    fig, ax = plt.subplots(figsize=(_FIG_W_IN, _FIG_H_IN), dpi=_FIG_DPI)
    fig.set_dpi(_FIG_DPI)
    _style_axes(ax)
    return fig, ax


def _title(view: _DistView, kind: str, *, include_run: bool = False) -> str:
    loc_bits = f"region={view.rname}, loc_id={view.loc}"
    if include_run and view.run_label:
        loc_bits = f"{view.run_label}; {loc_bits}"
    return (
        f"{view.fom} {kind} of {view.kernel_name}.{view.kernel_class} "
        f"on {view.host_label} ({loc_bits})"
    )


def _pct_bars_to_values(bars: Sequence[float], ref: float) -> List[float]:
    return [float(ref) * (1.0 + float(p) / 100.0) for p in bars]


def _value_xlim(
    values: np.ndarray,
    bars: Optional[Sequence[float]] = None,
    ref: Optional[float] = None,
) -> tuple:
    lo = float(np.min(values))
    hi = float(np.max(values))
    span = hi - lo
    pad = 1.0 if span <= 0.0 else max(1.0, 0.05 * span)
    x0 = lo - pad
    x1 = hi + pad
    if bars is not None and ref is not None:
        bar_vals = _pct_bars_to_values(bars, ref)
        x0 = min(x0, min(bar_vals))
        x1 = max(x1, max(bar_vals))
    if x1 <= x0:
        x1 = x0 + 2.0
    return (x0, x1)


def _data_xlim(arrays: Sequence[np.ndarray]) -> tuple:
    lo = min(float(np.min(a)) for a in arrays)
    hi = max(float(np.max(a)) for a in arrays)
    return _value_xlim(np.asarray([lo, hi], dtype=float))


def _resolve_axis_limit(
    user: Optional[tuple],
    auto: Optional[tuple],
) -> Optional[tuple]:
    if user is None:
        return auto
    if len(user) != 2:
        raise ValueError("axis limit must be a (min, max) tuple")
    lo = float(user[0])
    hi = float(user[1])
    if lo > hi:
        raise ValueError(
            f"axis limit min must be <= max, got ({lo}, {hi})"
        )
    return (lo, hi)


def _paint_gwr_guides(
    ax,
    bars: Sequence[float],
    cmap,
    n_bins: int,
    xlim,
    *,
    zero: float = 0.0,
) -> None:
    edges = np.concatenate(
        [[-np.inf], np.asarray(bars, dtype=float), [np.inf]]
    )
    x0, x1 = xlim
    for i in range(n_bins):
        lo = edges[i]
        hi = edges[i + 1]
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
    ax.axvline(float(zero), color="0.2", linestyle="-", linewidth=0.7, zorder=2)


def _attach_percent_top_axis(ax, xlim, ref: float, bars: Sequence[float]):
    if ref == 0.0:
        raise ValueError("ref=0; cannot form percent axis")
    vmin, vmax = xlim
    ax2 = ax.twiny()
    ax2.set_xlim((vmin / ref - 1.0) * 100.0, (vmax / ref - 1.0) * 100.0)
    ax2.set_xticks(list(bars))
    ax2.set_xticklabels([f"{p:g}%" for p in bars])
    ax2.set_xlabel("percent vs ref")
    _style_axes(ax2)
    return ax2


def _tab10_color(i: int):
    return plt.get_cmap("tab10")(i % 10)


def _group_views_by_region_loc(
    views: Sequence[_DistView],
) -> List[List[_DistView]]:
    groups: OrderedDict[tuple, List[_DistView]] = OrderedDict()
    for view in views:
        groups.setdefault((view.rid, view.loc), []).append(view)
    return list(groups.values())


def _collect_dist_views_for_root(
    data_root: str,
    kernel_name: str,
    kernel_class: str,
    *,
    hosts: Sequence[str],
    region_ids: Iterable[int],
    loc_ids: Iterable[int],
    ranks: Iterable[int],
    fom: str,
    xrange: Union[int, list, tuple],
    ref_section_base: bool,
    ref_key: str,
    ref_base: str,
    run_label: str,
) -> List[_DistView]:
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
            finite_ratio = ratio[np.isfinite(ratio)]
            if finite_ratio.size == 0:
                raise ValueError(
                    f"no finite value/ref samples for region_id={rid}, "
                    f"loc_id={loc} after filtering"
                )
            value_flat = draw[np.isfinite(draw)]
            if value_flat.size == 0:
                raise ValueError(
                    f"no finite {fom} samples for region_id={rid}, "
                    f"loc_id={loc} after filtering"
                )
            stat = draw if ref_section_base else mat
            finite_stat = stat[np.isfinite(stat)]
            mapping_ref = _compute_ref(finite_stat, ref_key)
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
                    n_samples=int(value_flat.size),
                    ratio_flat=np.asarray(finite_ratio, dtype=float),
                    value_flat=np.asarray(value_flat, dtype=float),
                    ref=float(mapping_ref),
                    run_label=run_label,
                )
            )

    if not views:
        raise ValueError(
            "no (region_id, loc_id) pairs to plot after filtering "
            f"(region_ids={selected_rids}, loc_ids={selected_locs}, "
            f"ranks={selected_ranks})"
        )
    return views


def _collect_dist_views(
    data_root: Union[str, Sequence[str]],
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

    roots = normalize_data_roots(data_root)
    multi_root = len(roots) > 1
    views: List[_DistView] = []
    for root in roots:
        run_label = os.path.basename(os.path.normpath(root))
        views.extend(
            _collect_dist_views_for_root(
                root,
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
                run_label=run_label,
            )
        )
    return views, list(bars), multi_root


def draw_histogram_npb_mpi(
    data_root: Union[str, Sequence[str]] = "",
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

    One figure per (region_id, loc_id) per data_root. Selected ranks are
    pooled into a single count histogram. X = percent-vs-ref bins from
    ``bars``; Y = count. Bar color is the same discrete green-white-red
    map as the heatmap.

    Parameters
    ----------
    data_root, kernel_name, kernel_class
        Path pieces for ``{data_root}/{kernel_name}.{kernel_class}/``.
        ``data_root`` may be a path string or a sequence of paths; each
        root is plotted independently. Remaining arguments are keyword-only.
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
    views, bars_used, multi_root = _collect_dist_views(
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
            ax.set_title(_title(view, "histogram", include_run=multi_root))
            fig.tight_layout()
            figures.append(fig)

    _display_notebook_figures(figures)
    return figures


def _draw_pdf_one(
    ax,
    view: _DistView,
    xlim: tuple,
    *,
    color,
    fill: bool,
    label: Optional[str],
) -> None:
    n_pdf_bins = int(np.clip(round(np.sqrt(view.value_flat.size)), 16, 80))
    dens, pdf_edges = np.histogram(
        view.value_flat, bins=n_pdf_bins, range=xlim, density=True
    )
    if fill:
        ax.stairs(
            dens,
            pdf_edges,
            fill=True,
            color=color,
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
            label=label,
        )
    else:
        ax.stairs(
            dens,
            pdf_edges,
            fill=True,
            color=color,
            alpha=0.18,
            linewidth=0.0,
            zorder=3,
        )
        ax.stairs(
            dens,
            pdf_edges,
            fill=False,
            color=color,
            linewidth=1.1,
            zorder=4,
            label=label,
        )


def draw_pdf_npb_mpi(
    data_root: Union[str, Sequence[str]] = "",
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
    overlap: bool = False,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
) -> List[Figure]:
    """
    Draw a pooled PDF of FOM values for NPB-MPI TacVar CSVs.

    One figure per (region_id, loc_id) when ``overlap`` is False (and per
    data_root). Bottom x = the FOM column; top x = percent-vs-ref ``bars``
    mapped through a single pooled ref. Selected ranks are pooled into one
    density histogram. No per-rank series.

    Parameters
    ----------
    data_root, kernel_name, kernel_class
        Path pieces for ``{data_root}/{kernel_name}.{kernel_class}/``.
        ``data_root`` may be a path string or a sequence of paths.
        Remaining arguments are keyword-only.
    hosts
        ``['all']`` / empty / None for all hosts; else short_host prefixes.
    region_ids, loc_ids
        Ids to plot; empty means all present (intersection that exists in data).
    ranks
        MPI ranks to include; empty means all ranks present. Values from
        selected ranks are flattened into one sample set (no per-rank series).
    fom
        Numeric CSV column (default ``elapsed_ns``).
    xrange
        Columns to include. ``-1`` (default) is all samples. A ``list`` of ints
        selects those indices. A ``tuple (a, b)`` of ints is the half-open
        slice ``[a, b)``; of floats in ``[0, 1]`` it is
        ``[int(n*a), int(n*b))``. Integer ``(0, 1)`` is only column 0.
    ref_section_base
        If True, compute the mapping ref on the ``xrange`` window; if False,
        on the full rank-by-sample matrix. Display still uses ``xrange``.
    ref_key
        ``max``, ``min``, ``median``, ``average``, or ``pxx`` percentile.
        Used for the percent top axis (pooled ``all`` aggregation).
    ref_base
        Sample set for per-cell ``value/ref`` used by the histogram path.
        PDF curves are raw FOM values; the twin percent axis always uses
        one pooled ref (``ref_key`` on the ``ref_section_base`` window).
    bars
        Strictly increasing percent-deviation cut points. Used as top-axis
        ticks and GWR background bands. Ignored when ``overlap`` is True.
    overlap
        If False (default), one figure per data_root. If True, overlay all
        listed roots on one axes per (region_id, loc_id) and ignore
        ``bars``.
    xlim
        Optional ``(min, max)`` for the main FOM x-axis (``min <= x <= max``).
        ``None`` (default) uses the auto range from data (and ``bars`` when
        ``overlap`` is False). Also used as the PDF histogram ``range``.
    ylim
        Optional ``(min, max)`` for the main density y-axis
        (``min <= y <= max``). ``None`` (default) leaves matplotlib autoscale.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    views, bars_used, multi_root = _collect_dist_views(
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

    figures: List[Figure] = []
    with plt.rc_context(_ACADEMIC_RC):
        if overlap:
            for group in _group_views_by_region_loc(views):
                xlim_used = _resolve_axis_limit(
                    xlim, _data_xlim([v.value_flat for v in group])
                )
                ylim_used = _resolve_axis_limit(ylim, None)
                fig, ax = _new_figure()
                for i, view in enumerate(group):
                    _draw_pdf_one(
                        ax,
                        view,
                        xlim_used,
                        color=_tab10_color(i),
                        fill=False,
                        label=view.run_label,
                    )
                ax.set_xlim(xlim_used)
                if ylim_used is not None:
                    ax.set_ylim(ylim_used)
                ax.set_xlabel(group[0].fom)
                ax.set_ylabel("probability density")
                ax.set_title(_title(group[0], "pdf", include_run=False))
                ax.legend()
                fig.tight_layout()
                figures.append(fig)
        else:
            edges, white_idx = _ratio_edges_from_bars(bars_used)
            n_bins = len(edges) - 1
            cmap = _green_white_red_cmap(n_bins, white_idx)
            for view in views:
                xlim_used = _resolve_axis_limit(
                    xlim, _value_xlim(view.value_flat, bars_used, view.ref)
                )
                ylim_used = _resolve_axis_limit(ylim, None)
                fig, ax = _new_figure()
                value_bars = _pct_bars_to_values(bars_used, view.ref)
                _paint_gwr_guides(
                    ax, value_bars, cmap, n_bins, xlim_used, zero=view.ref
                )
                _draw_pdf_one(
                    ax,
                    view,
                    xlim_used,
                    color="0.25",
                    fill=True,
                    label=None,
                )
                ax.set_xlim(xlim_used)
                if ylim_used is not None:
                    ax.set_ylim(ylim_used)
                ax.set_xlabel(view.fom)
                ax.set_ylabel("probability density")
                ax.set_title(_title(view, "pdf", include_run=multi_root))
                _attach_percent_top_axis(ax, xlim_used, view.ref, bars_used)
                fig.tight_layout()
                figures.append(fig)

    _display_notebook_figures(figures)
    return figures


def draw_cdf_npb_mpi(
    data_root: Union[str, Sequence[str]] = "",
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
    overlap: bool = False,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
) -> List[Figure]:
    """
    Draw a pooled empirical CDF of FOM values for NPB-MPI TacVar CSVs.

    One figure per (region_id, loc_id) when ``overlap`` is False (and per
    data_root). Bottom x = the FOM column; top x = percent-vs-ref ``bars``
    mapped through a single pooled ref. Selected ranks are pooled into one
    ECDF. A dotted line marks cumulative probability 0.5. No per-rank series.

    Parameters
    ----------
    data_root, kernel_name, kernel_class
        Path pieces for ``{data_root}/{kernel_name}.{kernel_class}/``.
        ``data_root`` may be a path string or a sequence of paths.
        Remaining arguments are keyword-only.
    hosts
        ``['all']`` / empty / None for all hosts; else short_host prefixes.
    region_ids, loc_ids
        Ids to plot; empty means all present (intersection that exists in data).
    ranks
        MPI ranks to include; empty means all ranks present. Values from
        selected ranks are flattened into one sample set (no per-rank series).
    fom
        Numeric CSV column (default ``elapsed_ns``).
    xrange
        Columns to include. ``-1`` (default) is all samples. A ``list`` of ints
        selects those indices. A ``tuple (a, b)`` of ints is the half-open
        slice ``[a, b)``; of floats in ``[0, 1]`` it is
        ``[int(n*a), int(n*b))``. Integer ``(0, 1)`` is only column 0.
    ref_section_base
        If True, compute the mapping ref on the ``xrange`` window; if False,
        on the full rank-by-sample matrix. Display still uses ``xrange``.
    ref_key
        ``max``, ``min``, ``median``, ``average``, or ``pxx`` percentile.
        Used for the percent top axis (pooled ``all`` aggregation).
    ref_base
        Sample set for per-cell ``value/ref`` used by the histogram path.
        CDF curves are raw FOM values; the twin percent axis always uses
        one pooled ref (``ref_key`` on the ``ref_section_base`` window).
    bars
        Strictly increasing percent-deviation cut points. Used as top-axis
        ticks and GWR background bands. Ignored when ``overlap`` is True.
    overlap
        If False (default), one figure per data_root. If True, overlay all
        listed roots on one axes per (region_id, loc_id) and ignore
        ``bars``.
    xlim
        Optional ``(min, max)`` for the main FOM x-axis (``min <= x <= max``).
        ``None`` (default) uses the auto range from data (and ``bars`` when
        ``overlap`` is False).
    ylim
        Optional ``(min, max)`` for the main cumulative-probability y-axis
        (``min <= y <= max``). ``None`` (default) uses ``(0.0, 1.02)``.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    views, bars_used, multi_root = _collect_dist_views(
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

    figures: List[Figure] = []
    with plt.rc_context(_ACADEMIC_RC):
        if overlap:
            for group in _group_views_by_region_loc(views):
                xlim_used = _resolve_axis_limit(
                    xlim, _data_xlim([v.value_flat for v in group])
                )
                ylim_used = _resolve_axis_limit(ylim, (0.0, 1.02))
                fig, ax = _new_figure()
                for i, view in enumerate(group):
                    xs = np.sort(view.value_flat)
                    ys = np.arange(1, xs.size + 1, dtype=float) / float(xs.size)
                    ax.step(
                        xs,
                        ys,
                        where="post",
                        color=_tab10_color(i),
                        linewidth=1.2,
                        zorder=4,
                        label=view.run_label,
                    )
                ax.axhline(0.5, color="0.5", linestyle=":", linewidth=0.8, zorder=3)
                ax.set_xlim(xlim_used)
                ax.set_ylim(ylim_used)
                ax.set_xlabel(group[0].fom)
                ax.set_ylabel("cumulative probability")
                ax.set_title(_title(group[0], "cdf", include_run=False))
                ax.legend()
                fig.tight_layout()
                figures.append(fig)
        else:
            edges, white_idx = _ratio_edges_from_bars(bars_used)
            n_bins = len(edges) - 1
            cmap = _green_white_red_cmap(n_bins, white_idx)
            for view in views:
                xlim_used = _resolve_axis_limit(
                    xlim, _value_xlim(view.value_flat, bars_used, view.ref)
                )
                ylim_used = _resolve_axis_limit(ylim, (0.0, 1.02))
                xs = np.sort(view.value_flat)
                ys = np.arange(1, xs.size + 1, dtype=float) / float(xs.size)
                fig, ax = _new_figure()
                value_bars = _pct_bars_to_values(bars_used, view.ref)
                _paint_gwr_guides(
                    ax, value_bars, cmap, n_bins, xlim_used, zero=view.ref
                )
                ax.step(xs, ys, where="post", color="black", linewidth=1.2, zorder=4)
                ax.axhline(0.5, color="0.5", linestyle=":", linewidth=0.8, zorder=3)
                ax.set_xlim(xlim_used)
                ax.set_ylim(ylim_used)
                ax.set_xlabel(view.fom)
                ax.set_ylabel("cumulative probability")
                ax.set_title(_title(view, "cdf", include_run=multi_root))
                _attach_percent_top_axis(ax, xlim_used, view.ref, bars_used)
                fig.tight_layout()
                figures.append(fig)

    _display_notebook_figures(figures)
    return figures
