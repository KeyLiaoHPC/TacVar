"""Generate campaign directories and guarded pretest/fulltest shell scripts."""
from __future__ import annotations

import json
import os
import shlex
import time
from dataclasses import replace
from pathlib import Path

from .spec import (
    CampaignSpec,
    digest_spec,
    ensure_kernel_np,
    read_spec,
    validate_spec,
    write_spec,
)

PKG_DIR = Path(__file__).resolve().parent
RUNNERS = PKG_DIR / "runners"
DEFAULT_WORKSPACE = Path(__file__).resolve().parents[3]
DEFAULT_SUITE = DEFAULT_WORKSPACE / "suites" / "NPB3.4.4" / "NPB3.4-MPI"


def _shell_quote_list(items: list[str] | tuple[str, ...]) -> str:
    return " ".join(shlex.quote(x) for x in items)


def _bash_array(name: str, items: list[str] | tuple[str, ...]) -> str:
    inner = " ".join(shlex.quote(x) for x in items)
    return f"{name}=({inner})"


def _render(template: str, mapping: dict[str, str]) -> str:
    out = template
    for key, val in mapping.items():
        out = out.replace(f"@{key}@", val)
    if "@" in out and any(f"@{k}@" in template for k in mapping):
        # leftover placeholders that look like ours
        leftover = [tok for tok in out.split() if tok.startswith("@") and tok.endswith("@")]
        if leftover:
            raise RuntimeError(f"unsubstituted placeholders: {leftover[:8]}")
    return out


def default_results_dir(spec: CampaignSpec, base: Path | None = None) -> Path:
    ts = time.strftime("%Y%m%dT%H%M%S")
    host = spec.hostname.replace("/", "_")
    root = base or Path(spec.suite_root or DEFAULT_SUITE)
    return root / f"results_{host}_class{spec.class_}_{ts}"


def materialize_spec(
    spec: CampaignSpec,
    results_dir: Path | None = None,
) -> tuple[Path, CampaignSpec, str]:
    spec = ensure_kernel_np(spec)
    errs = validate_spec(spec)
    if errs:
        raise ValueError("invalid campaign spec:\n- " + "\n- ".join(errs))
    if not spec.workspace:
        spec = replace(spec, workspace=str(DEFAULT_WORKSPACE))
    if not spec.suite_root:
        spec = replace(spec, suite_root=str(DEFAULT_SUITE))
    results = Path(results_dir) if results_dir else default_results_dir(spec)
    results.mkdir(parents=True, exist_ok=True)
    (results / "preflight").mkdir(exist_ok=True)
    (results / "runs").mkdir(exist_ok=True)
    (results / "summary" / "figures").mkdir(parents=True, exist_ok=True)
    write_spec(spec, results)
    dig = digest_spec(spec)
    return results, spec, dig


def _profiles_json(spec: CampaignSpec) -> str:
    return json.dumps([p.to_dict() for p in spec.counter_profiles], separators=(",", ":"))


def _kernel_np_json(spec: CampaignSpec) -> str:
    return json.dumps(spec.kernel_np, separators=(",", ":"), sort_keys=True)


def _common_map(spec: CampaignSpec, results_dir: Path, dig: str) -> dict[str, str]:
    return {
        "HOSTNAME": spec.hostname,
        "ARCH": spec.arch,
        "CORES": str(spec.cores),
        "CLASS": spec.class_,
        "MPI_HOME": spec.mpi_home,
        "PAPI_HOME": spec.papi_home,
        "WORKSPACE": spec.workspace,
        "SUITE_ROOT": spec.suite_root,
        "RESULTS_DIR": str(results_dir),
        "SPEC_DIGEST": dig,
        "BIND_ARGS": _shell_quote_list(spec.bind_args),
        "BIND_ARRAY": _bash_array("BIND_ARGS", spec.bind_args),
        "KERNELS_ARRAY": _bash_array("KERNELS", spec.kernels),
        "TIMERS_ARRAY": _bash_array("TIMERS", spec.timers),
        "PROFILES_JSON": _profiles_json(spec),
        "KERNEL_NP_JSON": _kernel_np_json(spec),
        "PRETEST_NP": str(spec.pretest_np),
        "PRETEST_KERNEL": spec.pretest_kernel,
        "PRETEST_CLASS": spec.pretest_class,
        "TACVAR_NSTP": str(spec.tacvar_nstp),
        "PER_STEP": "1" if spec.per_step else "0",
        "NPB_TIMER_FLAG": "1" if spec.npb_timer_flag else "0",
        "LOCK_PATH": f"/tmp/tacvar-npb-mpi-timer-{os.environ.get('USER', 'user')}.lock",
    }


def generate_pretest(
    spec: CampaignSpec,
    results_dir: Path | None = None,
) -> Path:
    results, spec, dig = materialize_spec(spec, results_dir)
    tpl = (RUNNERS / "pretest.sh.in").read_text()
    script = _render(tpl, _common_map(spec, results, dig))
    out = results / "00_pretest.sh"
    out.write_text(script)
    out.chmod(out.stat().st_mode | 0o111)
    analyze = results / "02_analyze.sh"
    analyze.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f'RESULTS={shlex.quote(str(results))}\n'
        f'SCRIPTS={shlex.quote(str(DEFAULT_WORKSPACE / "utils" / "scripts"))}\n'
        'cd "$SCRIPTS"\n'
        'if [[ -x .venv/bin/python ]]; then PY=.venv/bin/python; else PY=python3; fi\n'
        '"$PY" analyze_npb_timer_campaign.py "$RESULTS" "$@"\n'
    )
    analyze.chmod(analyze.stat().st_mode | 0o111)
    return out


def _preflight_ok(results_dir: Path, dig: str) -> None:
    status = results_dir / "preflight" / "STATUS"
    sha = results_dir / "preflight" / "spec.sha256"
    if not status.exists() or status.read_text().strip() != "PASS":
        raise RuntimeError(
            f"pretest not PASS (see {status}). Refuse to generate fulltest."
        )
    if not sha.exists():
        raise RuntimeError(f"missing {sha}")
    got = sha.read_text().strip()
    if got != dig:
        raise RuntimeError(
            f"pretest digest mismatch: preflight={got} campaign={dig}. "
            "Update the specification and rerun pretest."
        )


def generate_fulltest(results_dir: Path) -> Path:
    results_dir = Path(results_dir).resolve()
    spec = read_spec(results_dir / "campaign.json")
    dig = digest_spec(spec)
    campaign_sha = (results_dir / "campaign.sha256").read_text().strip()
    if campaign_sha != dig:
        raise RuntimeError("campaign.sha256 does not match campaign.json")
    _preflight_ok(results_dir, dig)
    tpl = (RUNNERS / "fulltest.sh.in").read_text()
    script = _render(tpl, _common_map(spec, results_dir, dig))
    out = results_dir / "01_fulltest.sh"
    out.write_text(script)
    out.chmod(out.stat().st_mode | 0o111)
    return out


def load_spec_file(path: Path) -> CampaignSpec:
    return read_spec(Path(path))
