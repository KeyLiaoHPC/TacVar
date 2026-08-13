"""Discrete green-white-red heatmaps for NPB-MPI TacVar measurements."""

from __future__ import annotations

import os
import re
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgb
from matplotlib.figure import Figure

from .io import (
    load_measurement_frame,
    load_timer_info_table,
    normalize_data_roots,
    resolve_kernel_dir,
    select_ids,
    short_hosts_from_frame,
)

_PXX_RE = re.compile(r"^p(\d{1,3})$", re.IGNORECASE)

_DEFAULT_BARS = [-15, -5, -1, 1, 5, 15]

_DEEP_GREEN = to_rgb("#00441b")
_LIGHT_GREEN = to_rgb("#a1d99b")
_WHITE = to_rgb("#ffffff")
_LIGHT_RED = to_rgb("#fcbba1")
_DEEP_RED = to_rgb("#cb181d")

_ACADEMIC_RC = {
    "font.family": "serif",
    "font.serif": [
        "Cambria",
        "Liberation Serif",
        "DejaVu Serif",
        "Times New Roman",
        "serif",
    ],
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 100,
}


def _compute_ref(values: np.ndarray, ref_key: str) -> float:
    key = (ref_key or "").strip().lower()
    if values.size == 0:
        raise ValueError("cannot compute ref: empty sample set")
    if key == "max":
        ref = float(np.max(values))
    elif key == "min":
        ref = float(np.min(values))
    elif key == "median":
        ref = float(np.median(values))
    elif key in ("average", "mean"):
        ref = float(np.mean(values))
    else:
        m = _PXX_RE.match(key)
        if not m:
            raise ValueError(
                f"unsupported ref_key={ref_key!r}; "
                "use max, min, median, average, or pxx (e.g. p50, p90)"
            )
        pct = int(m.group(1))
        if pct < 0 or pct > 100:
            raise ValueError(f"percentile in ref_key must be in [0, 100], got p{pct}")
        ref = float(np.percentile(values, pct))
    if ref == 0.0:
        raise ValueError(f"ref_key={ref_key!r} produced ref=0; cannot form value/ref")
    return ref


def _is_bool(x) -> bool:
    return isinstance(x, (bool, np.bool_))


def _is_int_like(x) -> bool:
    return isinstance(x, (int, np.integer)) and not _is_bool(x)


def _is_float_like(x) -> bool:
    return isinstance(x, (float, np.floating)) and not _is_bool(x)


def _resolve_xrange(n: int, xrange) -> np.ndarray:
    """
    Map xrange to column indices in ``[0, n)``.

    ``-1``: all columns. ``list``: those integer indices (Python negatives
    allowed). ``tuple (a, b)``: integer slice ``[a, b)``, or fractional
    ``[int(n*a), int(n*b))`` when both values are floats in ``[0, 1]``.
    """
    if n <= 0:
        raise ValueError("cannot resolve xrange: no sample columns")
    if xrange == -1:
        return np.arange(n, dtype=int)
    if isinstance(xrange, list):
        if not xrange:
            raise ValueError("xrange list is empty")
        idx = []
        for v in xrange:
            if not _is_int_like(v):
                raise ValueError(
                    f"xrange list entries must be integers, got {v!r}"
                )
            iv = int(v)
            if iv < 0:
                iv = n + iv
            if iv < 0 or iv >= n:
                raise ValueError(
                    f"xrange index {v} out of range for n={n} columns"
                )
            idx.append(iv)
        return np.asarray(idx, dtype=int)
    if isinstance(xrange, tuple):
        if len(xrange) != 2:
            raise ValueError(
                f"xrange tuple must have length 2, got {xrange!r}"
            )
        a, b = xrange
        if _is_int_like(a) and _is_int_like(b):
            ia, ib = int(a), int(b)
            if ia > ib:
                raise ValueError(
                    f"xrange tuple requires a <= b, got ({ia}, {ib})"
                )
            selected = np.arange(n, dtype=int)[ia:ib]
        elif _is_float_like(a) and _is_float_like(b):
            fa, fb = float(a), float(b)
            if not (0.0 <= fa <= fb <= 1.0):
                raise ValueError(
                    "xrange float tuple requires 0 <= a <= b <= 1, "
                    f"got ({a}, {b})"
                )
            selected = np.arange(int(n * fa), int(n * fb), dtype=int)
        else:
            raise ValueError(
                "xrange tuple must be two ints [a, b) or two floats in "
                f"[0, 1], got {xrange!r}"
            )
        if selected.size == 0:
            raise ValueError(
                f"xrange={xrange!r} selected no columns from n={n}"
            )
        return selected
    raise ValueError(
        "xrange must be -1, a list of indices, or a (a, b) tuple, "
        f"got {type(xrange).__name__}"
    )


def _ratio_from_refs(
    mat: np.ndarray,
    col_idx: np.ndarray,
    ref_key: str,
    ref_base: str,
    ref_section_base: bool,
) -> np.ndarray:
    """
    ``value / ref`` for drawn columns ``col_idx``.

    Stats come from the xrange window when ``ref_section_base`` is True,
    otherwise from the full matrix. ``ref_base`` is ``all``, ``x``, or ``y``.
    """
    base = (ref_base or "all").strip().lower()
    draw = mat[:, col_idx]
    stat = draw if ref_section_base else mat
    ratio = np.full(draw.shape, np.nan, dtype=float)

    if base == "all":
        finite = stat[np.isfinite(stat)]
        ref = _compute_ref(finite, ref_key)
        np.divide(draw, ref, out=ratio)
        return ratio

    if base == "x":
        for j, orig_j in enumerate(col_idx):
            col = mat[:, orig_j]
            finite = col[np.isfinite(col)]
            if finite.size == 0:
                continue
            ref = _compute_ref(finite, ref_key)
            ratio[:, j] = draw[:, j] / ref
        return ratio

    if base == "y":
        for i in range(mat.shape[0]):
            row = stat[i, :]
            finite = row[np.isfinite(row)]
            if finite.size == 0:
                continue
            ref = _compute_ref(finite, ref_key)
            ratio[i, :] = draw[i, :] / ref
        return ratio

    raise ValueError(
        f"unsupported ref_base={ref_base!r}; use all, x, or y"
    )


def _ratio_edges_from_bars(bars: Sequence[float]) -> Tuple[np.ndarray, int]:
    """
    Build edges for len(bars)+1 bins.

    bars are percent deviations. Convert to ratio edges (1 + p/100), then
    add +/- inf. Returns (edges of length n_bins+1, index of the bin that
    contains ratio 1.0).
    """
    if not bars:
        raise ValueError("bars must be a non-empty sequence of percent cut points")

    ratios = np.asarray([1.0 + float(p) / 100.0 for p in bars], dtype=float)
    if np.any(np.diff(ratios) <= 0):
        raise ValueError(f"bars must be strictly increasing, got {list(bars)}")

    edges = np.concatenate([[-np.inf], ratios, [np.inf]])
    cut_finite = edges[1:-1]
    white_idx = int(np.digitize(1.0, cut_finite, right=True))
    return edges, white_idx


def _lerp(c0, c1, t: float):
    t = float(np.clip(t, 0.0, 1.0))
    return tuple(c0[i] + (c1[i] - c0[i]) * t for i in range(3))


def _green_white_red_cmap(n_bins: int, white_idx: int) -> ListedColormap:
    """
    Discrete colormap: more negative vs ref → deeper green;
    0%-containing bin → white; more positive → deeper red.
    """
    if n_bins < 3:
        raise ValueError(f"need at least 3 bins, got {n_bins}")
    if white_idx < 0 or white_idx >= n_bins:
        raise ValueError(f"white_idx={white_idx} out of range for n_bins={n_bins}")

    colors = []
    n_neg = white_idx
    n_pos = n_bins - white_idx - 1
    for i in range(n_bins):
        if i < white_idx:
            # leftmost = deepest green; near white = light green
            if n_neg == 1:
                t = 0.0
            else:
                t = i / (n_neg - 1)
            colors.append(_lerp(_DEEP_GREEN, _LIGHT_GREEN, t))
        elif i == white_idx:
            colors.append(_WHITE)
        else:
            # near white = light red; rightmost = deepest red
            j = i - white_idx - 1
            if n_pos == 1:
                t = 1.0
            else:
                t = j / (n_pos - 1)
            colors.append(_lerp(_LIGHT_RED, _DEEP_RED, t))
    return ListedColormap(colors, name="tacvar_gwr")


def _build_rank_matrix(
    sub: pd.DataFrame, fom: str
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build (n_ranks, n_samples) matrix of fom values.
    Rows are sorted ranks; columns are sample indices in CSV order per rank.
    Missing cells are NaN.
    """
    if fom not in sub.columns:
        raise ValueError(
            f"fom column {fom!r} not in CSV; columns={list(sub.columns)}"
        )
    ranks = sorted({int(r) for r in sub["rank"].unique()})
    series_by_rank = {}
    max_len = 0
    for rank in ranks:
        vals = sub.loc[sub["rank"] == rank, fom].to_numpy(dtype=float)
        series_by_rank[rank] = vals
        max_len = max(max_len, len(vals))
    if max_len == 0:
        raise ValueError("no samples for selected region/loc")
    mat = np.full((len(ranks), max_len), np.nan, dtype=float)
    for i, rank in enumerate(ranks):
        vals = series_by_rank[rank]
        mat[i, : len(vals)] = vals
    return mat, np.asarray(ranks, dtype=int)


def _format_edge_label(edge: float) -> str:
    if np.isneginf(edge):
        return "-inf"
    if np.isposinf(edge):
        return "+inf"
    pct = (edge - 1.0) * 100.0
    if abs(pct) < 1e-9:
        return "0%"
    return f"{pct:g}%"


def _format_host_label(hosts: Sequence[str]) -> str:
    if not hosts:
        return "all"
    if len(hosts) == 1:
        return hosts[0]
    return ",".join(hosts)


def _display_notebook_figures(figures: Sequence[Figure]) -> None:
    """Render each figure in Jupyter/VS Code.

    A cell whose last expression is a list of Figure objects only prints
    ``[<Figure ...>]``; the inline backend often skips a second show on
    rerun. Explicit ``display`` makes each rerun paint the heatmaps.
    """
    try:
        ip = get_ipython()  # type: ignore[name-defined]
    except NameError:
        return
    if ip is None:
        return
    try:
        from IPython.display import display
    except ImportError:
        return
    for fig in figures:
        display(fig)


def _heatmap_figures_for_root(
    data_root: str,
    kernel_name: str,
    kernel_class: str,
    *,
    hosts: Sequence[str],
    region_ids: Iterable[int],
    loc_ids: Iterable[int],
    fom: str,
    xrange: Union[int, list, tuple],
    ref_section_base: bool,
    ref_key: str,
    ref_base: str,
    cmap: ListedColormap,
    cut_points: np.ndarray,
    norm: BoundaryNorm,
    n_bins: int,
    run_label: str,
) -> List[Figure]:
    """Draw heatmaps for one data_root. Does not display figures."""
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
    host_label = _format_host_label(short_hosts_from_frame(df))

    figures: List[Figure] = []
    for rid in selected_rids:
        rname = region_names.get(rid, f"r{rid}")
        for loc in selected_locs:
            sub = df[(df["region_id"] == rid) & (df["loc_id"] == loc)]
            if sub.empty:
                continue
            mat, ranks = _build_rank_matrix(sub, fom)
            col_idx = _resolve_xrange(int(mat.shape[1]), xrange)
            draw = mat[:, col_idx]
            n_ranks = len(ranks)
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

            bin_idx = np.full(ratio.shape, np.nan, dtype=float)
            finite_mask = np.isfinite(ratio)
            bin_idx[finite_mask] = np.digitize(
                ratio[finite_mask], cut_points, right=True
            ).astype(float)

            width_in = max(8.0, n_pts * 0.12 + 3.0)
            height_in = 10.24  # 1024 px at dpi=100
            fig, ax = plt.subplots(figsize=(width_in, height_in), dpi=100)
            fig.set_dpi(100)

            cmap_with_bad = cmap.copy()
            cmap_with_bad.set_bad(color="#e0e0e0")
            im = ax.imshow(
                bin_idx,
                aspect="auto",
                interpolation="nearest",
                cmap=cmap_with_bad,
                norm=norm,
                origin="upper",
            )
            for spine in ax.spines.values():
                spine.set_color("black")
                spine.set_linewidth(0.8)

            ax.set_xlabel("sample index (CSV order)")
            ax.set_ylabel("MPI rank")
            ax.set_yticks(np.arange(len(ranks)))
            ax.set_yticklabels([str(r) for r in ranks])
            if n_pts <= 40:
                xtick_pos = np.arange(n_pts)
            else:
                step = max(1, n_pts // 20)
                xtick_pos = np.arange(0, n_pts, step)
            ax.set_xticks(xtick_pos)
            ax.set_xticklabels([str(int(col_idx[i])) for i in xtick_pos])

            loc_bits = f"region={rname}, loc_id={loc}"
            if run_label:
                loc_bits = f"{run_label}; {loc_bits}"
            title = (
                f"{fom} heatmap of {kernel_name}.{kernel_class} "
                f"on {host_label} ({loc_bits})"
            )
            ax.set_title(title)

            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_ticks(np.arange(n_bins - 1) + 0.5)
            cbar.set_ticklabels([_format_edge_label(e) for e in cut_points])
            cbar.set_label("percent vs ref (green < 0 < red)")

            fig.tight_layout()
            figures.append(fig)

    if not figures:
        raise ValueError(
            "no (region_id, loc_id) pairs to plot after filtering "
            f"(region_ids={selected_rids}, loc_ids={selected_locs})"
        )
    return figures


def draw_heatmap_npb_mpi(
    data_root: Union[str, Sequence[str]] = "",
    kernel_name: str = "",
    kernel_class: str = "",
    *,
    hosts: Optional[Sequence[str]] = None,
    region_ids: Optional[Iterable[int]] = None,
    loc_ids: Optional[Iterable[int]] = None,
    fom: str = "elapsed_ns",
    xrange: Union[int, list, tuple] = -1,
    ref_section_base: bool = True,
    ref_key: str = "median",
    ref_base: str = "all",
    bars: Optional[Sequence[float]] = None,
) -> List[Figure]:
    """
    Draw discrete green-white-red heatmaps for NPB-MPI TacVar CSVs.

    One figure per (region_id, loc_id) per data_root. X = sample index
    (CSV order), Y = MPI rank. Cell color is the discrete bin of value/ref.

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
    fom
        Numeric CSV column (default ``elapsed_ns``).
    xrange
        Columns to draw. ``-1`` (default) is all samples. A ``list`` of ints
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
        ``ref_section_base``).
    bars
        Strictly increasing percent-deviation cut points. Number of bins is
        ``len(bars)+1`` after wrapping with ``+/- inf``. The bin that
        contains ratio 1.0 is white.

    Returns
    -------
    list of matplotlib.figure.Figure
    """
    if hosts is None:
        hosts = ["all"]
    if region_ids is None:
        region_ids = range(0)
    if loc_ids is None:
        loc_ids = range(0)
    if bars is None:
        bars = list(_DEFAULT_BARS)

    roots = normalize_data_roots(data_root)
    multi_root = len(roots) > 1
    edges, white_idx = _ratio_edges_from_bars(bars)
    n_bins = len(edges) - 1
    cmap = _green_white_red_cmap(n_bins, white_idx)
    cut_points = edges[1:-1]
    bin_edges_cb = np.arange(n_bins + 1) - 0.5
    norm = BoundaryNorm(bin_edges_cb, ncolors=n_bins)

    figures: List[Figure] = []
    with plt.rc_context(_ACADEMIC_RC):
        for root in roots:
            run_label = (
                os.path.basename(os.path.normpath(root)) if multi_root else ""
            )
            figures.extend(
                _heatmap_figures_for_root(
                    root,
                    kernel_name,
                    kernel_class,
                    hosts=hosts,
                    region_ids=region_ids,
                    loc_ids=loc_ids,
                    fom=fom,
                    xrange=xrange,
                    ref_section_base=ref_section_base,
                    ref_key=ref_key,
                    ref_base=ref_base,
                    cmap=cmap,
                    cut_points=cut_points,
                    norm=norm,
                    n_bins=n_bins,
                    run_label=run_label,
                )
            )

    _display_notebook_figures(figures)
    return figures
