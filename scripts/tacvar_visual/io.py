"""Load NPB-MPI TacVar measurement CSVs for visualization."""

from __future__ import annotations

import glob
import os
import re
from typing import Iterable, List, Optional, Sequence

import pandas as pd

_CSV_NAME_RE = re.compile(
    r"^(?P<host>.+)_r(?P<rank>\d{4})_t(?P<thread>\d{4})_p(?P<pid>\d+)\.csv$"
)


def resolve_kernel_dir(
    data_root: str,
    kernel_name: str,
    kernel_class: str,
) -> str:
    """Resolve ``{data_root}/{kernel_name}.{kernel_class}/``."""
    if not data_root or not kernel_name or not kernel_class:
        raise ValueError(
            "data_root, kernel_name, and kernel_class are required "
            f"(got data_root={data_root!r}, kernel_name={kernel_name!r}, "
            f"kernel_class={kernel_class!r})"
        )

    base = os.path.join(data_root, f"{kernel_name}.{kernel_class}")
    if not os.path.isdir(base):
        raise ValueError(f"kernel directory not found: {base!r}")
    return base


def list_rank_csvs(kernel_dir: str, hosts: Optional[Sequence[str]] = None) -> List[str]:
    """List per-rank CSV paths under kernel_dir, optionally filtered by short_host."""
    pattern = os.path.join(kernel_dir, "*_r*_t*_p*.csv")
    paths = sorted(glob.glob(pattern))
    selected: List[str] = []
    host_filter = _normalize_hosts(hosts)
    for path in paths:
        name = os.path.basename(path)
        m = _CSV_NAME_RE.match(name)
        if not m:
            continue
        if host_filter is not None and m.group("host") not in host_filter:
            continue
        selected.append(path)
    if not selected:
        raise ValueError(
            f"no rank CSV files matching '*_r*_t*_p*.csv' under {kernel_dir!r}"
            + (f" for hosts={list(hosts)!r}" if host_filter is not None else "")
        )
    return selected


def _normalize_hosts(hosts: Optional[Sequence[str]]) -> Optional[set]:
    if hosts is None:
        return None
    host_list = list(hosts)
    if not host_list:
        return None
    non_str = [h for h in host_list if not isinstance(h, str)]
    if non_str:
        raise ValueError(
            "hosts must be short_host strings such as ['all'] or ['c920bn1'], "
            f"not region/loc ids; got {non_str!r}. Pass hosts=, region_ids=, "
            "and loc_ids= as keywords, then re-import tacvar_visual."
        )
    if len(host_list) == 1 and host_list[0] in ("all", ""):
        return None
    if all(h in ("all", "") for h in host_list):
        return None
    return {h for h in host_list if h and h != "all"}


def short_hosts_from_frame(df: pd.DataFrame) -> List[str]:
    """Unique short_host prefixes from ``_source_file``, sorted."""
    if "_source_file" not in df.columns:
        return []
    hosts = []
    for name in df["_source_file"].astype(str).unique():
        m = _CSV_NAME_RE.match(name)
        if m:
            hosts.append(m.group("host"))
    return sorted(set(hosts))


def load_timer_info_table(kernel_dir: str) -> pd.DataFrame:
    """Load timer_info.csv as a DataFrame; empty if missing."""
    path = os.path.join(kernel_dir, "timer_info.csv")
    if not os.path.isfile(path):
        return pd.DataFrame()
    return pd.read_csv(path)


def load_timer_info(kernel_dir: str) -> dict:
    """Load region_id -> name from timer_info.csv; empty dict if missing."""
    df = load_timer_info_table(kernel_dir)
    if df.empty or "region_id" not in df.columns or "name" not in df.columns:
        return {}
    return {int(r): str(n) for r, n in zip(df["region_id"], df["name"])}


def load_measurement_frame(
    kernel_dir: str,
    hosts: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Concatenate selected rank CSVs into one DataFrame."""
    paths = list_rank_csvs(kernel_dir, hosts)
    frames = []
    for path in paths:
        df = pd.read_csv(path)
        df["_source_file"] = os.path.basename(path)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def select_ids(df: pd.DataFrame, column: str, wanted: Iterable[int]) -> List[int]:
    """Return sorted ids from ``column``; empty ``wanted`` means all present."""
    if column not in df.columns:
        raise ValueError(f"column {column!r} not in frame")
    available = sorted({int(x) for x in df[column].unique()})
    wanted_list = list(wanted)
    if not wanted_list:
        return available
    wanted_set = {int(x) for x in wanted_list}
    selected = [r for r in available if r in wanted_set]
    if not selected:
        raise ValueError(
            f"none of requested {column}={sorted(wanted_set)} found; "
            f"available={available}"
        )
    return selected
