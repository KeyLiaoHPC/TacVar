"""Histogram and inverse-CDF plotting for NPB-MPI per-step timings."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np

from .load import Campaign, collect_step_by_timer, load_campaign
from .plot_style import LINE_STYLES, MARKERS, apply_academic_style
from .statistics import choose_bin_edges, histogram_fraction, percentile, quantiles_ms, summarize


def _as_campaign(obj: Campaign | str | Path) -> Campaign:
    if isinstance(obj, Campaign):
        return obj
    return load_campaign(obj)


def _filter_cells(
    campaign: Campaign,
    kernels: Iterable[str] | None,
    timers: Iterable[str] | None,
    profiles: Iterable[str] | None,
) -> tuple[list[str], list[str], list[str]]:
    ks = sorted({c.kernel for c in campaign.cells})
    ts = sorted({c.timer for c in campaign.cells})
    ps = sorted({c.profile for c in campaign.cells})
    if kernels is not None:
        want = {k.lower() for k in kernels}
        ks = [k for k in ks if k in want]
    if timers is not None:
        want = set(timers)
        ts = [t for t in ts if t in want]
    if profiles is not None:
        want = set(profiles)
        ps = [p for p in ps if p in want]
    return ks, ts, ps


def _savefig(fig: plt.Figure, base: Path, formats: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    base.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        p = base.with_suffix(f".{fmt}")
        kw: dict[str, Any] = {}
        if fmt == "png":
            kw["dpi"] = 300
        fig.savefig(p, bbox_inches="tight", **kw)
        paths.append(p)
    plt.close(fig)
    return paths


def plot_step_histograms(
    campaign: Campaign | str | Path,
    *,
    out_dir: Path | None = None,
    formats: Sequence[str] = ("png", "pdf", "eps"),
    kernels: Iterable[str] | None = None,
    timers: Iterable[str] | None = None,
    profiles: Iterable[str] | None = None,
    allow_font_fallback: bool = False,
    nbins: int = 24,
) -> dict[str, Any]:
    camp = _as_campaign(campaign)
    style_meta = apply_academic_style(allow_font_fallback=allow_font_fallback)
    fig_root = out_dir or (camp.results_dir / "summary" / "figures")
    fig_root.mkdir(parents=True, exist_ok=True)
    ks, ts, ps = _filter_cells(camp, kernels, timers, profiles)
    payload: dict[str, Any] = {
        "style": style_meta,
        "kernels": ks,
        "timers": ts,
        "profiles": ps,
        "figures": [],
        "series": {},
    }
    paths: list[Path] = []
    for kernel in ks:
        for profile in ps:
            by_timer = collect_step_by_timer(camp, kernel, profile, ts)
            if not by_timer:
                continue
            all_ms = [v / 1e6 for vals in by_timer.values() for v in vals]
            edges, scale = choose_bin_edges(all_ms, nbins=nbins)
            fig, ax = plt.subplots(figsize=(5.2, 3.6))
            for i, timer in enumerate(sorted(by_timer)):
                vals_ms = [v / 1e6 for v in by_timer[timer]]
                frac = histogram_fraction(vals_ms, edges)
                centers = [
                    0.5 * (edges[j] + edges[j + 1]) for j in range(len(edges) - 1)
                ]
                ax.step(
                    centers,
                    frac,
                    where="mid",
                    linestyle=LINE_STYLES[i % len(LINE_STYLES)],
                    color="black",
                    linewidth=1.0,
                    label=timer,
                )
                payload["series"][f"{kernel}/{profile}/{timer}"] = {
                    "n": len(vals_ms),
                    "unit": "ms",
                    "edges": edges,
                    "fraction": frac,
                    "summary": summarize(by_timer[timer]),
                }
            ax.set_xlabel("Per-step runtime (ms)")
            ax.set_ylabel("Fraction of samples")
            ax.set_title(f"{kernel.upper()} / {profile} per-step histogram")
            if scale == "log":
                ax.set_xscale("log")
            ax.legend(frameon=False, loc="best")
            ax.text(
                0.0,
                -0.22,
                f"Source: {camp.results_dir.name} | bins={scale}",
                transform=ax.transAxes,
                fontsize=7,
                va="top",
            )
            base = fig_root / f"hist_{kernel}_{profile}"
            fig_paths = _savefig(fig, base, formats)
            paths.extend(fig_paths)
            payload["figures"].append(
                {
                    "kernel": kernel,
                    "profile": profile,
                    "scale": scale,
                    "edges": edges,
                    "files": [str(p.relative_to(camp.results_dir)) for p in fig_paths],
                }
            )
    json_path = camp.results_dir / "summary" / "step_histograms.json"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    payload["json"] = str(json_path)
    payload["paths"] = [str(p) for p in paths]
    return payload


def plot_step_icdf(
    campaign: Campaign | str | Path,
    *,
    out_dir: Path | None = None,
    formats: Sequence[str] = ("png", "pdf", "eps"),
    kernels: Iterable[str] | None = None,
    timers: Iterable[str] | None = None,
    profiles: Iterable[str] | None = None,
    allow_font_fallback: bool = False,
    full_tail: bool = False,
    pct_grid: Sequence[float] | None = None,
) -> dict[str, Any]:
    camp = _as_campaign(campaign)
    style_meta = apply_academic_style(allow_font_fallback=allow_font_fallback)
    fig_root = out_dir or (camp.results_dir / "summary" / "figures")
    fig_root.mkdir(parents=True, exist_ok=True)
    ks, ts, ps = _filter_cells(camp, kernels, timers, profiles)
    if pct_grid is None:
        pct_grid = list(np.linspace(0, 100, 201))
    payload: dict[str, Any] = {
        "style": style_meta,
        "kernels": ks,
        "timers": ts,
        "profiles": ps,
        "full_tail": full_tail,
        "figures": [],
        "series": {},
    }
    paths: list[Path] = []
    for kernel in ks:
        for profile in ps:
            by_timer = collect_step_by_timer(camp, kernel, profile, ts)
            if not by_timer:
                continue
            fig, ax = plt.subplots(figsize=(5.2, 3.6))
            xmax = 100.0 if full_tail else 99.0
            for i, timer in enumerate(sorted(by_timer)):
                q = quantiles_ms(by_timer[timer], pct_grid)
                payload["series"][f"{kernel}/{profile}/{timer}"] = {
                    "n": len(by_timer[timer]),
                    "percentiles": list(pct_grid),
                    "runtime_ms": q,
                    "summary": summarize(by_timer[timer]),
                }
                xs = [p for p in pct_grid if p <= xmax]
                ys = [q[j] for j, p in enumerate(pct_grid) if p <= xmax]
                ax.plot(
                    xs,
                    ys,
                    linestyle=LINE_STYLES[i % len(LINE_STYLES)],
                    marker=MARKERS[i % len(MARKERS)],
                    markevery=max(len(xs) // 12, 1),
                    markersize=3.5,
                    color="black",
                    linewidth=1.0,
                    label=timer,
                )
            ax.set_xlabel("Percentile (%)")
            ax.set_ylabel("Per-step runtime (ms)")
            title_tail = "full tail" if full_tail else "p99-focused"
            ax.set_title(f"{kernel.upper()} / {profile} inverse CDF ({title_tail})")
            ax.set_xlim(0, xmax)
            ax.legend(frameon=False, loc="best")
            ax.text(
                0.0,
                -0.22,
                f"Source: {camp.results_dir.name}",
                transform=ax.transAxes,
                fontsize=7,
                va="top",
            )
            base = fig_root / f"icdf_{kernel}_{profile}"
            fig_paths = _savefig(fig, base, formats)
            paths.extend(fig_paths)
            payload["figures"].append(
                {
                    "kernel": kernel,
                    "profile": profile,
                    "files": [str(p.relative_to(camp.results_dir)) for p in fig_paths],
                }
            )
    json_path = camp.results_dir / "summary" / "step_icdf.json"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    payload["json"] = str(json_path)
    payload["paths"] = [str(p) for p in paths]
    return payload
