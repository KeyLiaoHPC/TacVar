"""Numeric helpers for campaign summaries and plots."""
from __future__ import annotations

import math
import statistics
from typing import Sequence


def percentile(sorted_vals: Sequence[float], p: float) -> float:
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    k = (len(sorted_vals) - 1) * (p / 100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return float(sorted_vals[int(k)])
    return float(sorted_vals[f] * (c - k) + sorted_vals[c] * (k - f))


def summarize(vals: Sequence[float]) -> dict[str, float | int]:
    if not vals:
        return {
            "n": 0,
            "min": float("nan"),
            "p50": float("nan"),
            "mean": float("nan"),
            "p90": float("nan"),
            "p99": float("nan"),
            "max": float("nan"),
            "cv": float("nan"),
        }
    s = sorted(float(v) for v in vals)
    mean = statistics.fmean(s)
    return {
        "n": len(s),
        "min": s[0],
        "p50": percentile(s, 50),
        "mean": mean,
        "p90": percentile(s, 90),
        "p99": percentile(s, 99),
        "max": s[-1],
        "cv": (statistics.pstdev(s) / mean) if mean else float("nan"),
    }


def quantiles_ms(vals_ns: Sequence[float], percentiles: Sequence[float]) -> list[float]:
    if not vals_ns:
        return [float("nan")] * len(percentiles)
    s = sorted(float(v) / 1e6 for v in vals_ns)
    return [percentile(s, p) for p in percentiles]


def choose_bin_edges(all_ms: Sequence[float], nbins: int = 24) -> tuple[list[float], str]:
    if not all_ms:
        return [0.0, 1.0], "linear"
    s = sorted(float(v) for v in all_ms)
    lo = percentile(s, 1)
    hi = percentile(s, 99)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        lo, hi = s[0], s[-1] if s[-1] > s[0] else s[0] + 1.0
    use_log = lo > 0 and (hi / lo) > 20
    if use_log:
        llo, lhi = math.log10(lo), math.log10(hi)
        edges = [10 ** (llo + (lhi - llo) * i / nbins) for i in range(nbins + 1)]
        return edges, "log"
    edges = [lo + (hi - lo) * i / nbins for i in range(nbins + 1)]
    return edges, "linear"


def histogram_fraction(vals_ms: Sequence[float], edges: Sequence[float]) -> list[float]:
    nbins = len(edges) - 1
    counts = [0] * nbins
    if not vals_ms or nbins <= 0:
        return [0.0] * max(nbins, 0)
    for v in vals_ms:
        if v <= edges[0]:
            counts[0] += 1
            continue
        if v >= edges[-1]:
            counts[-1] += 1
            continue
        idx = nbins - 1
        for i in range(nbins):
            if edges[i] <= v < edges[i + 1] or (i == nbins - 1 and v <= edges[i + 1]):
                idx = i
                break
        counts[idx] += 1
    n = float(len(vals_ms))
    return [c / n for c in counts]
