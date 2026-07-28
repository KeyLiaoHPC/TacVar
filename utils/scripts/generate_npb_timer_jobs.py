#!/usr/bin/env python3
"""Generate NPB-MPI timer campaign pretest/fulltest scripts."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow running without install: python3 utils/scripts/generate_npb_timer_jobs.py
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from npb_timer_campaign.generate import generate_fulltest, generate_pretest
from npb_timer_campaign.spec import (
    CampaignSpec,
    CounterProfile,
    ensure_kernel_np,
    nearest_square_le,
    read_spec,
    validate_spec,
    write_spec,
)


def _parse_profiles(raw: str) -> tuple[CounterProfile, ...]:
    data = json.loads(raw)
    return tuple(CounterProfile.from_dict(p) for p in data)


def cmd_init(args: argparse.Namespace) -> int:
    profiles = _parse_profiles(args.profiles_json)
    kernels = tuple(k.upper() for k in args.kernels.split(","))
    timers = tuple(t.strip() for t in args.timers.split(",") if t.strip())
    knp = {}
    if args.kernel_np_json:
        knp = {k.upper(): int(v) for k, v in json.loads(args.kernel_np_json).items()}
    for k in kernels:
        if k in ("BT", "SP") and k not in knp:
            if args.square_np:
                knp[k] = int(args.square_np)
            else:
                knp[k] = nearest_square_le(args.cores)
    spec = CampaignSpec(
        hostname=args.hostname,
        cores=args.cores,
        kernels=kernels,
        timers=timers,
        counter_profiles=profiles,
        class_=args.class_,
        mpi_home=args.mpi_home,
        papi_home=args.papi_home,
        bind_args=tuple(args.bind_args.split()) if args.bind_args else ("--map-by", "core", "--bind-to", "core"),
        kernel_np=knp,
        workspace=args.workspace or "",
        suite_root=args.suite_root or "",
        arch=args.arch,
        pretest_np=args.pretest_np,
        tacvar_nstp=args.tacvar_nstp,
    )
    spec = ensure_kernel_np(spec)
    errs = validate_spec(spec)
    if errs:
        print("invalid spec:\n- " + "\n- ".join(errs), file=sys.stderr)
        return 2
    results = Path(args.results_dir) if args.results_dir else None
    out = generate_pretest(spec, results)
    print(out.parent)
    print(out)
    return 0


def cmd_pretest(args: argparse.Namespace) -> int:
    spec_path = Path(args.spec).resolve()
    spec = read_spec(spec_path)
    if args.results_dir:
        results = Path(args.results_dir)
    else:
        # If campaign.json already lives in a results tree, reuse that directory.
        results = spec_path.parent
    out = generate_pretest(spec, results)
    print(out.parent)
    print(out)
    return 0


def cmd_fulltest(args: argparse.Namespace) -> int:
    out = generate_fulltest(Path(args.results_dir))
    print(out)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser("init", help="Create results dir + campaign.json + 00_pretest.sh")
    p_init.add_argument("--hostname", required=True)
    p_init.add_argument("--cores", type=int, required=True)
    p_init.add_argument("--kernels", required=True, help="Comma-separated, e.g. BT,CG,EP")
    p_init.add_argument("--timers", required=True, help="Comma-separated timer names")
    p_init.add_argument("--profiles-json", required=True, help="JSON array of counter profiles")
    p_init.add_argument("--class", dest="class_", required=True)
    p_init.add_argument("--mpi-home", required=True)
    p_init.add_argument("--papi-home", required=True)
    p_init.add_argument("--bind-args", default="--map-by core --bind-to core")
    p_init.add_argument("--kernel-np-json", default="")
    p_init.add_argument("--square-np", type=int, default=0, help="NP for BT/SP when cores not square")
    p_init.add_argument("--workspace", default="")
    p_init.add_argument("--suite-root", default="")
    p_init.add_argument("--results-dir", default="")
    p_init.add_argument("--arch", default="aarch64")
    p_init.add_argument("--pretest-np", type=int, default=4)
    p_init.add_argument("--tacvar-nstp", type=int, default=10)
    p_init.set_defaults(func=cmd_init)

    p_pre = sub.add_parser("pretest", help="Write 00_pretest.sh from an existing campaign.json")
    p_pre.add_argument("--spec", required=True)
    p_pre.add_argument("--results-dir", default="")
    p_pre.set_defaults(func=cmd_pretest)

    p_full = sub.add_parser("fulltest", help="Write 01_fulltest.sh after matching pretest PASS")
    p_full.add_argument("--results-dir", required=True)
    p_full.set_defaults(func=cmd_fulltest)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
