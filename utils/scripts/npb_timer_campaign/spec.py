"""Campaign specification: validation, canonical JSON, and SHA-256 digest."""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

VALID_KERNELS = ("BT", "CG", "EP", "FT", "IS", "LU", "MG", "SP")
VALID_CLASSES = ("S", "W", "A", "B", "C", "D", "E")
VALID_TIMERS_AARCH64 = (
    "native",
    "mpi_wtime",
    "clock_gettime",
    "papi_get_real_nsec",
    "cntvct_el0",
    "cntvct_el0_dmb",
)
VALID_TIMERS_X86 = (
    "native",
    "mpi_wtime",
    "clock_gettime",
    "papi_get_real_nsec",
    "rdtsc",
    "rdtscp",
    "rdtscp_lfence",
    "tsc_asym",
)
VALID_COUNTER_BACKENDS = ("none", "papi_read", "perf_event_open", "asm")
SQUARE_KERNELS = ("BT", "SP")


@dataclass(frozen=True)
class CounterProfile:
    name: str
    backend: str
    events: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "backend": self.backend,
            "events": list(self.events),
        }

    @staticmethod
    def from_dict(d: dict[str, Any]) -> "CounterProfile":
        events = tuple(d.get("events") or ())
        return CounterProfile(
            name=str(d["name"]),
            backend=str(d["backend"]),
            events=events,
        )


@dataclass(frozen=True)
class CampaignSpec:
    hostname: str
    cores: int
    kernels: tuple[str, ...]
    timers: tuple[str, ...]
    counter_profiles: tuple[CounterProfile, ...]
    class_: str
    mpi_home: str
    papi_home: str
    bind_args: tuple[str, ...] = ("--map-by", "core", "--bind-to", "core")
    kernel_np: dict[str, int] = field(default_factory=dict)
    workspace: str = ""
    suite_root: str = ""
    per_step: bool = True
    npb_timer_flag: bool = True
    pretest_np: int = 4
    pretest_kernel: str = "MG"
    pretest_class: str = "S"
    tacvar_nstp: int = 10
    arch: str = "aarch64"

    def to_dict(self) -> dict[str, Any]:
        return {
            "hostname": self.hostname,
            "cores": self.cores,
            "kernels": list(self.kernels),
            "timers": list(self.timers),
            "counter_profiles": [p.to_dict() for p in self.counter_profiles],
            "class": self.class_,
            "mpi_home": self.mpi_home,
            "papi_home": self.papi_home,
            "bind_args": list(self.bind_args),
            "kernel_np": dict(self.kernel_np),
            "workspace": self.workspace,
            "suite_root": self.suite_root,
            "per_step": self.per_step,
            "npb_timer_flag": self.npb_timer_flag,
            "pretest_np": self.pretest_np,
            "pretest_kernel": self.pretest_kernel,
            "pretest_class": self.pretest_class,
            "tacvar_nstp": self.tacvar_nstp,
            "arch": self.arch,
        }

    @staticmethod
    def from_dict(d: dict[str, Any]) -> "CampaignSpec":
        profiles = tuple(
            CounterProfile.from_dict(p) for p in (d.get("counter_profiles") or [])
        )
        kernels = tuple(k.upper() for k in (d.get("kernels") or []))
        return CampaignSpec(
            hostname=str(d["hostname"]),
            cores=int(d["cores"]),
            kernels=kernels,
            timers=tuple(d.get("timers") or []),
            counter_profiles=profiles,
            class_=str(d.get("class") or d.get("class_")),
            mpi_home=str(d["mpi_home"]),
            papi_home=str(d["papi_home"]),
            bind_args=tuple(d.get("bind_args") or ("--map-by", "core", "--bind-to", "core")),
            kernel_np={k.upper(): int(v) for k, v in (d.get("kernel_np") or {}).items()},
            workspace=str(d.get("workspace") or ""),
            suite_root=str(d.get("suite_root") or ""),
            per_step=bool(d.get("per_step", True)),
            npb_timer_flag=bool(d.get("npb_timer_flag", True)),
            pretest_np=int(d.get("pretest_np", 4)),
            pretest_kernel=str(d.get("pretest_kernel", "MG")).upper(),
            pretest_class=str(d.get("pretest_class", "S")),
            tacvar_nstp=int(d.get("tacvar_nstp", 10)),
            arch=str(d.get("arch", "aarch64")),
        )


def is_square(n: int) -> bool:
    if n <= 0:
        return False
    r = int(math.isqrt(n))
    return r * r == n


def nearest_square_le(n: int) -> int:
    if n <= 0:
        raise ValueError("cores must be positive")
    r = int(math.isqrt(n))
    return r * r


def validate_spec(spec: CampaignSpec) -> list[str]:
    errs: list[str] = []
    if not spec.hostname.strip():
        errs.append("hostname is required")
    if spec.cores <= 0:
        errs.append("cores must be a positive integer")
    if not spec.kernels:
        errs.append("kernels must be non-empty")
    for k in spec.kernels:
        if k not in VALID_KERNELS:
            errs.append(f"invalid kernel: {k}")
    if not spec.timers:
        errs.append("timers must be non-empty")
    allowed = VALID_TIMERS_AARCH64 if spec.arch == "aarch64" else VALID_TIMERS_X86
    for t in spec.timers:
        if t not in allowed:
            errs.append(f"timer {t} not valid for arch={spec.arch}")
    if not spec.counter_profiles:
        errs.append("counter_profiles must be non-empty")
    for p in spec.counter_profiles:
        if p.backend not in VALID_COUNTER_BACKENDS:
            errs.append(f"invalid counter backend: {p.backend}")
        if p.backend == "none" and p.events:
            errs.append(f"profile {p.name}: none backend requires empty events")
        if p.backend != "none" and not p.events:
            errs.append(f"profile {p.name}: non-none backend requires events")
        if p.backend != "none" and len(p.events) == 0:
            errs.append(f"profile {p.name}: events empty")
    if spec.class_ not in VALID_CLASSES:
        errs.append(f"invalid class: {spec.class_}")
    if not spec.mpi_home:
        errs.append("mpi_home is required")
    if not spec.papi_home:
        errs.append("papi_home is required")
    for k in spec.kernels:
        np = spec.kernel_np.get(k, spec.cores)
        if k in SQUARE_KERNELS and not is_square(np):
            errs.append(
                f"kernel {k} requires square process count; got {np}. "
                "Set kernel_np or omit the kernel."
            )
    return errs


def canonical_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def digest_spec(spec: CampaignSpec) -> str:
    return hashlib.sha256(canonical_json(spec.to_dict()).encode("utf-8")).hexdigest()


def write_spec(spec: CampaignSpec, results_dir: Path) -> tuple[Path, Path]:
    errs = validate_spec(spec)
    if errs:
        raise ValueError("invalid campaign spec:\n- " + "\n- ".join(errs))
    results_dir.mkdir(parents=True, exist_ok=True)
    json_path = results_dir / "campaign.json"
    sha_path = results_dir / "campaign.sha256"
    payload = spec.to_dict()
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    dig = digest_spec(spec)
    sha_path.write_text(dig + "\n")
    return json_path, sha_path


def read_spec(path: Path) -> CampaignSpec:
    data = json.loads(path.read_text())
    spec = CampaignSpec.from_dict(data)
    errs = validate_spec(spec)
    if errs:
        raise ValueError("invalid campaign.spec:\n- " + "\n- ".join(errs))
    return spec


def ensure_kernel_np(spec: CampaignSpec) -> CampaignSpec:
    """Fill missing kernel_np with cores (caller must already resolve square cases)."""
    knp = dict(spec.kernel_np)
    for k in spec.kernels:
        knp.setdefault(k, spec.cores)
    return CampaignSpec(
        hostname=spec.hostname,
        cores=spec.cores,
        kernels=spec.kernels,
        timers=spec.timers,
        counter_profiles=spec.counter_profiles,
        class_=spec.class_,
        mpi_home=spec.mpi_home,
        papi_home=spec.papi_home,
        bind_args=spec.bind_args,
        kernel_np=knp,
        workspace=spec.workspace,
        suite_root=spec.suite_root,
        per_step=spec.per_step,
        npb_timer_flag=spec.npb_timer_flag,
        pretest_np=spec.pretest_np,
        pretest_kernel=spec.pretest_kernel,
        pretest_class=spec.pretest_class,
        tacvar_nstp=spec.tacvar_nstp,
        arch=spec.arch,
    )
