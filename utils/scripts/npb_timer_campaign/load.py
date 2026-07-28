"""Discover and load NPB-MPI campaign result trees (campaign.json or legacy)."""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator

from .spec import CampaignSpec, CounterProfile, read_spec

NAS_TIME_RE = re.compile(r"Time in seconds\s*=\s*([0-9.]+)", re.I)
NAS_CLASS_RE = re.compile(r"Class\s*=\s*(\S+)", re.I)
NAS_PROCS_RE = re.compile(r"Total processes\s*=\s*(\d+)", re.I)
DETAIL_TIMER_RE = re.compile(
    r"^\s*timer\s+(\d+)\s*\(\s*([^)]+?)\s*\)\s*:?\s*"
    r"([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)\s*$",
    re.I,
)
CONF_RE = {
    "timer": re.compile(r"^TACVAR_TIMER=(.*)$", re.M),
    "backend": re.compile(r"^TACVAR_COUNTER_BACKEND=(.*)$", re.M),
    "count": re.compile(r"^TACVAR_COUNTER_COUNT=(.*)$", re.M),
    "names": re.compile(r"^TACVAR_COUNTER_NAMES=(.*)$", re.M),
    "nstp": re.compile(r"^TACVAR_NSTP=(.*)$", re.M),
    "per_step": re.compile(r"^TACVAR_ENABLE_PER_STEP_TIMING=(.*)$", re.M),
}


@dataclass
class RunCell:
    timer: str
    profile: str
    kernel: str
    path: Path
    status: str = ""
    backend: str = ""
    events: tuple[str, ...] = ()
    nstp: int = 0
    nas_time_s: float | None = None
    nas_class: str = ""
    nas_np: int | None = None
    data_dir: Path | None = None


@dataclass
class Campaign:
    results_dir: Path
    spec: CampaignSpec | None = None
    cells: list[RunCell] = field(default_factory=list)
    meta: dict[str, str] = field(default_factory=dict)
    events_file: str = ""
    legacy: bool = False

    @property
    def kernels(self) -> list[str]:
        return sorted({c.kernel for c in self.cells})

    @property
    def timers(self) -> list[str]:
        return sorted({c.timer for c in self.cells})

    @property
    def profiles(self) -> list[str]:
        return sorted({c.profile for c in self.cells})


def _parse_conf(text: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, rx in CONF_RE.items():
        m = rx.search(text)
        out[key] = m.group(1).strip() if m else ""
    names = out.get("names") or ""
    out["events"] = tuple(x for x in names.split(",") if x) if names else ()
    try:
        out["nstp"] = int(out.get("nstp") or 0)
    except ValueError:
        out["nstp"] = 0
    return out


def _find_data_dir(run_dir: Path) -> Path | None:
    dirs = sorted(run_dir.glob("data_????????T??????"))
    return dirs[-1] if dirs else None


def _parse_stdout(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    text = path.read_text(errors="replace")
    out: dict[str, Any] = {"successful": "SUCCESSFUL" in text}
    m = NAS_TIME_RE.search(text)
    out["nas_time_s"] = float(m.group(1)) if m else None
    m = NAS_CLASS_RE.search(text)
    out["nas_class"] = m.group(1) if m else ""
    m = NAS_PROCS_RE.search(text)
    out["nas_np"] = int(m.group(1)) if m else None
    details = []
    for line in text.splitlines():
        dm = DETAIL_TIMER_RE.match(line)
        if dm:
            details.append(
                {
                    "timer_id": int(dm.group(1)),
                    "timer_name": re.sub(r"\s+", "", dm.group(2).strip()),
                    "min_s": float(dm.group(3)),
                    "max_s": float(dm.group(4)),
                    "avg_s": float(dm.group(5)),
                }
            )
    out["detail"] = details
    return out


def _split_tag(tag: str) -> tuple[str, str]:
    # native_e0 / clock_gettime_e4 / cntvct_el0_dmb_e0
    m = re.match(r"^(.*)_(e\d+)$", tag)
    if not m:
        return tag, "unknown"
    return m.group(1), m.group(2)


def discover_legacy(results_dir: Path) -> Campaign:
    cells: list[RunCell] = []
    runs = results_dir / "runs"
    for tag_dir in sorted(runs.glob("*")) if runs.is_dir() else []:
        if not tag_dir.is_dir():
            continue
        timer, profile = _split_tag(tag_dir.name)
        for kern_dir in sorted(tag_dir.iterdir()):
            if not kern_dir.is_dir():
                continue
            conf_p = kern_dir / "tacvar.conf"
            conf = _parse_conf(conf_p.read_text()) if conf_p.exists() else {}
            st = (kern_dir / "STATUS").read_text().strip() if (kern_dir / "STATUS").exists() else ""
            so = _parse_stdout(kern_dir / "stdout.log")
            cells.append(
                RunCell(
                    timer=timer,
                    profile=profile,
                    kernel=kern_dir.name.lower(),
                    path=kern_dir,
                    status=st,
                    backend=str(conf.get("backend") or ""),
                    events=tuple(conf.get("events") or ()),
                    nstp=int(conf.get("nstp") or 0),
                    nas_time_s=so.get("nas_time_s"),
                    nas_class=str(so.get("nas_class") or ""),
                    nas_np=so.get("nas_np"),
                    data_dir=_find_data_dir(kern_dir),
                )
            )
    meta: dict[str, str] = {}
    meta_p = results_dir / "meta.env"
    if meta_p.exists():
        for line in meta_p.read_text().splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                meta[k] = v
    events = ""
    ef = results_dir / "events4.txt"
    if ef.exists():
        events = ef.read_text().strip()
    return Campaign(
        results_dir=results_dir,
        spec=None,
        cells=cells,
        meta=meta,
        events_file=events,
        legacy=True,
    )


def discover_campaign(results_dir: Path) -> Campaign:
    results_dir = Path(results_dir).resolve()
    spec_path = results_dir / "campaign.json"
    if spec_path.exists():
        spec = read_spec(spec_path)
        # Still discover cells from runs/ for completed data.
        camp = discover_legacy(results_dir)
        camp.spec = spec
        camp.legacy = False
        return camp
    return discover_legacy(results_dir)


def load_campaign(results_dir: str | Path) -> Campaign:
    return discover_campaign(Path(results_dir))


def iter_step_elapsed_ns(cell: RunCell) -> Iterator[float]:
    if not cell.data_dir:
        return
        yield  # pragma: no cover
    for csv_path in sorted(cell.data_dir.glob("npb-mpi_*.csv")):
        with csv_path.open(newline="") as f:
            for row in csv.DictReader(f):
                try:
                    if int(row.get("region_id", "-1")) == 1000:
                        yield float(row["elapsed_ns"])
                except (TypeError, ValueError, KeyError):
                    continue


def iter_total_elapsed_ns(cell: RunCell) -> Iterator[float]:
    if not cell.data_dir:
        return
        yield  # pragma: no cover
    for csv_path in sorted(cell.data_dir.glob("npb-mpi_*.csv")):
        with csv_path.open(newline="") as f:
            for row in csv.DictReader(f):
                try:
                    if int(row.get("region_id", "-1")) != 1000:
                        yield float(row["elapsed_ns"])
                except (TypeError, ValueError, KeyError):
                    continue


def collect_step_by_timer(
    campaign: Campaign,
    kernel: str,
    profile: str,
    timers: Iterable[str] | None = None,
) -> dict[str, list[float]]:
    want = set(timers) if timers is not None else None
    out: dict[str, list[float]] = {}
    for cell in campaign.cells:
        if cell.kernel != kernel or cell.profile != profile:
            continue
        if want is not None and cell.timer not in want:
            continue
        vals = list(iter_step_elapsed_ns(cell))
        if vals:
            out[cell.timer] = vals
    return out


def parse_detail_from_stdout(stdout_path: Path) -> list[dict[str, Any]]:
    return list(_parse_stdout(stdout_path).get("detail") or [])
