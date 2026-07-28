"""Summarize campaign runs and write durable REPORT.md."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from .load import (
    Campaign,
    RunCell,
    iter_step_elapsed_ns,
    iter_total_elapsed_ns,
    load_campaign,
    parse_detail_from_stdout,
)
from .spec import digest_spec
from .statistics import summarize


def _as_campaign(obj: Campaign | str | Path) -> Campaign:
    if isinstance(obj, Campaign):
        return obj
    return load_campaign(obj)


def _matrix_complete(camp: Campaign) -> tuple[bool, list[str]]:
    if camp.spec is None:
        # Legacy: all discovered cells should be PASS
        bad = [f"{c.timer}_{c.profile}/{c.kernel}:{c.status}" for c in camp.cells if c.status != "PASS"]
        return (not bad), bad
    missing: list[str] = []
    for timer in camp.spec.timers:
        for p in camp.spec.counter_profiles:
            for k in camp.spec.kernels:
                tag = f"{timer}_{p.name}"
                hit = [
                    c
                    for c in camp.cells
                    if c.timer == timer and c.profile == p.name and c.kernel == k.lower()
                ]
                if not hit:
                    missing.append(f"{tag}/{k}:MISSING")
                elif hit[0].status != "PASS":
                    missing.append(f"{tag}/{k}:{hit[0].status}")
    return (not missing), missing


def summarize_campaign(campaign: Campaign | str | Path) -> dict[str, Any]:
    camp = _as_campaign(campaign)
    summary_dir = camp.results_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    overall_rows: list[dict[str, Any]] = []
    total_rows: list[dict[str, Any]] = []
    step_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []

    for cell in camp.cells:
        overall_rows.append(
            {
                "timer": cell.timer,
                "profile": cell.profile,
                "kernel": cell.kernel,
                "status": cell.status,
                "nas_time_s": cell.nas_time_s if cell.nas_time_s is not None else "",
                "nas_class": cell.nas_class,
                "nas_np": cell.nas_np if cell.nas_np is not None else "",
                "backend": cell.backend,
                "events": ",".join(cell.events),
            }
        )
        tot = list(iter_total_elapsed_ns(cell))
        if tot:
            st = summarize(tot)
            total_rows.append(
                {
                    "timer": cell.timer,
                    "profile": cell.profile,
                    "kernel": cell.kernel,
                    **{k: st[k] for k in ("n", "min", "p50", "mean", "p90", "p99", "max", "cv")},
                }
            )
        step = list(iter_step_elapsed_ns(cell))
        if step:
            st = summarize(step)
            step_rows.append(
                {
                    "timer": cell.timer,
                    "profile": cell.profile,
                    "kernel": cell.kernel,
                    **{k: st[k] for k in ("n", "min", "p50", "mean", "p90", "p99", "max", "cv")},
                }
            )
        stdout = cell.path / "stdout.log"
        for d in parse_detail_from_stdout(stdout):
            detail_rows.append(
                {
                    "timer": cell.timer,
                    "profile": cell.profile,
                    "kernel": cell.kernel,
                    **d,
                }
            )

    def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
        if not rows:
            path.write_text("")
            return
        keys = list(rows[0].keys())
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)

    write_csv(summary_dir / "overall_runtime.csv", overall_rows)
    write_csv(summary_dir / "total_dist_stats.csv", total_rows)
    write_csv(summary_dir / "step_dist_stats.csv", step_rows)
    write_csv(summary_dir / "detail_from_stdout.csv", detail_rows)

    complete, missing = _matrix_complete(camp)
    analysis = {
        "kernels": camp.kernels,
        "timers": camp.timers,
        "profiles": camp.profiles,
        "n_cells": len(camp.cells),
        "n_pass": sum(1 for c in camp.cells if c.status == "PASS"),
        "complete": complete,
        "missing_or_fail": missing,
        "legacy": camp.legacy,
        "meta": camp.meta,
        "events_file": camp.events_file,
    }
    # Data-driven conclusions: stable kernels = low CV across timers for e0 step p50
    conclusions = _derive_conclusions(camp, step_rows)
    analysis["conclusions"] = conclusions
    (summary_dir / "step_dist_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
    return analysis


def _derive_conclusions(camp: Campaign, step_rows: list[dict[str, Any]]) -> list[str]:
    out: list[str] = []
    by_kp: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in step_rows:
        by_kp[(r["kernel"], r["profile"])].append(r)
    for (kernel, profile), rows in sorted(by_kp.items()):
        p50s = [float(r["p50"]) for r in rows if r.get("p50") == r.get("p50")]
        if len(p50s) < 2:
            continue
        mean = sum(p50s) / len(p50s)
        if mean <= 0:
            continue
        spread = (max(p50s) - min(p50s)) / mean
        cvs = [float(r["cv"]) for r in rows if r.get("cv") == r.get("cv")]
        med_cv = sorted(cvs)[len(cvs) // 2] if cvs else float("nan")
        if spread < 0.05:
            out.append(
                f"{kernel.upper()}/{profile}: timer p50s agree within {spread:.1%} "
                f"(timer-stable for this matrix)."
            )
        elif spread > 0.25:
            out.append(
                f"{kernel.upper()}/{profile}: timer p50s disagree by {spread:.1%} "
                f"(timer-sensitive)."
            )
        if med_cv == med_cv and med_cv > 0.3:
            out.append(
                f"{kernel.upper()}/{profile}: heavy-tailed per-step CV (median CV={med_cv:.2f})."
            )
    # e0 vs nonzero profile overall NAS time deltas when both present
    profiles = camp.profiles
    if len(profiles) >= 2:
        base = profiles[0]
        for other in profiles[1:]:
            deltas = []
            for k in camp.kernels:
                t0 = [
                    c.nas_time_s
                    for c in camp.cells
                    if c.kernel == k and c.profile == base and c.nas_time_s
                ]
                t1 = [
                    c.nas_time_s
                    for c in camp.cells
                    if c.kernel == k and c.profile == other and c.nas_time_s
                ]
                if t0 and t1:
                    # average across timers
                    a0 = sum(t0) / len(t0)
                    a1 = sum(t1) / len(t1)
                    if a0:
                        deltas.append((k, (a1 - a0) / a0))
            if deltas:
                mean_d = sum(d for _, d in deltas) / len(deltas)
                out.append(
                    f"Mean NAS-time change {base}->{other}: {mean_d:+.1%} "
                    f"across {len(deltas)} kernels (not a uniform tax)."
                )
    if not out:
        out.append("No strong cross-timer pattern detected for the discovered matrix.")
    return out


def write_report(
    campaign: Campaign | str | Path,
    summary: dict[str, Any] | None = None,
    hist_paths: dict[str, Any] | None = None,
    icdf_paths: dict[str, Any] | None = None,
    *,
    partial: bool = False,
) -> Path:
    camp = _as_campaign(campaign)
    if summary is None:
        summary = summarize_campaign(camp)
    report = camp.results_dir / "REPORT.md"
    dig = ""
    spec_block = "_Legacy result tree (no campaign.json). Discovery inferred the matrix._\n"
    if camp.spec is not None:
        dig = digest_spec(camp.spec)
        spec_block = (
            "```json\n"
            + json.dumps(camp.spec.to_dict(), indent=2, sort_keys=True)
            + "\n```\n"
            + f"\n**Digest:** `{dig}`\n"
        )

    status_line = "PARTIAL" if partial else ("COMPLETE" if summary.get("complete") else "INCOMPLETE")
    lines: list[str] = []
    lines.append("# NPB-MPI Timer Comparison Report\n")
    lines.append(f"**Status:** {status_line}\n")
    lines.append(f"**Results:** `{camp.results_dir}`\n")
    lines.append("\n## Campaign specification\n\n")
    lines.append(spec_block)
    lines.append("\n## Host / toolchain\n\n")
    if camp.meta:
        lines.append("| Key | Value |\n|---|---|\n")
        for k, v in sorted(camp.meta.items()):
            lines.append(f"| `{k}` | `{v}` |\n")
    else:
        lines.append("_No meta.env present._\n")
    if camp.events_file:
        lines.append(f"\n**Events file:** `{camp.events_file}`\n")

    lines.append("\n## Matrix completion\n\n")
    lines.append(
        f"- Cells discovered: **{summary.get('n_cells', 0)}** "
        f"(PASS={summary.get('n_pass', 0)})\n"
    )
    lines.append(f"- Kernels: {', '.join(summary.get('kernels') or [])}\n")
    lines.append(f"- Timers: {', '.join(summary.get('timers') or [])}\n")
    lines.append(f"- Profiles: {', '.join(summary.get('profiles') or [])}\n")
    miss = summary.get("missing_or_fail") or []
    if miss:
        lines.append("\nMissing/FAIL cells:\n\n")
        for m in miss:
            lines.append(f"- `{m}`\n")

    lines.append("\n## Conclusions\n\n")
    for c in summary.get("conclusions") or []:
        lines.append(f"- {c}\n")

    lines.append("\n## Summary tables\n\n")
    for name in (
        "overall_runtime.csv",
        "total_dist_stats.csv",
        "detail_from_stdout.csv",
        "step_dist_stats.csv",
    ):
        p = camp.results_dir / "summary" / name
        if p.exists() and p.stat().st_size > 0:
            lines.append(f"- [`summary/{name}`](summary/{name})\n")

    lines.append("\n## Figures\n\n")
    fig_dir = camp.results_dir / "summary" / "figures"
    if fig_dir.is_dir():
        pngs = sorted(fig_dir.glob("*.png"))
        for png in pngs:
            rel = png.relative_to(camp.results_dir).as_posix()
            stem = png.stem
            lines.append(f"### {stem}\n\n")
            lines.append(f"![{stem}]({rel})\n\n")
            for ext in ("pdf", "eps"):
                alt = png.with_suffix(f".{ext}")
                if alt.exists():
                    lines.append(f"- [{ext.upper()}]({alt.relative_to(camp.results_dir).as_posix()})\n")
            lines.append("\n")
    else:
        lines.append("_No figures yet._\n")

    if hist_paths and hist_paths.get("style"):
        st = hist_paths["style"]
        lines.append("\n## Plot style\n\n")
        lines.append(
            f"- Requested font: `{st.get('requested_font')}`\n"
            f"- Resolved font: `{st.get('resolved_font')}`\n"
            f"- Fallback used: `{st.get('font_fallback')}`\n"
        )

    lines.append("\n## Limitations\n\n")
    if partial:
        lines.append("- Report generated with `--allow-partial`; treat as diagnostic.\n")
    if camp.legacy:
        lines.append("- Legacy tree without `campaign.json`; matrix inferred from `runs/`.\n")
    lines.append("- Raw data under `runs/<timer>_<profile>/<kernel>/data_*/`.\n")
    lines.append("\n")

    report.write_text("".join(lines))
    return report


def analyze_campaign(
    results_dir: str | Path,
    *,
    formats: tuple[str, ...] = ("png", "pdf", "eps"),
    allow_partial: bool = False,
    allow_font_fallback: bool = False,
    kernels: Iterable[str] | None = None,
    timers: Iterable[str] | None = None,
    profiles: Iterable[str] | None = None,
    full_tail: bool = False,
) -> dict[str, Any]:
    from .plots import plot_step_histograms, plot_step_icdf

    camp = load_campaign(results_dir)
    summary = summarize_campaign(camp)
    complete = bool(summary.get("complete"))
    if not complete and not allow_partial:
        raise RuntimeError(
            "Campaign matrix incomplete; refuse analysis. "
            "Pass allow_partial=True / --allow-partial for diagnostics."
        )
    hist = plot_step_histograms(
        camp,
        formats=formats,
        kernels=kernels,
        timers=timers,
        profiles=profiles,
        allow_font_fallback=allow_font_fallback,
    )
    icdf = plot_step_icdf(
        camp,
        formats=formats,
        kernels=kernels,
        timers=timers,
        profiles=profiles,
        allow_font_fallback=allow_font_fallback,
        full_tail=full_tail,
    )
    report = write_report(
        camp,
        summary,
        hist,
        icdf,
        partial=allow_partial and not complete,
    )
    return {
        "campaign": camp,
        "summary": summary,
        "histograms": hist,
        "icdf": icdf,
        "report": report,
    }
