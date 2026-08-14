#!/usr/bin/env python3
"""
@file gen_config.py
@brief Generate tacvar_generated_config.h/.mk from tacvar.conf (build-time).

Detects target arch via: compiler -dM -E -
Validates timer/counter/consumer combinations; emits only selected backends.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

CONSUMERS = {
    "npb-mpi": "TACVAR_CONSUMER_NPB_MPI",
    "npb-omp": "TACVAR_CONSUMER_NPB_OMP",
    "lmbench": "TACVAR_CONSUMER_LMBENCH",
}

TIMERS_ALL = [
    "native",
    "gettimeofday",
    "omp_get_wtime",
    "mpi_wtime",
    "clock_gettime",
    "papi_get_real_nsec",
    "rdtsc",
    "rdtscp",
    "rdtscp_lfence",
    "tsc_asym",
    "cntvct_el0",
    "cntvct_el0_dmb",
]

TIMERS_X86 = {
    "rdtsc", "rdtscp", "rdtscp_lfence", "tsc_asym",
}
TIMERS_ARM = {"cntvct_el0", "cntvct_el0_dmb"}
TIMERS_TICK = TIMERS_X86 | TIMERS_ARM

COUNTERS = ["none", "perf_event_open", "papi_read", "asm"]

TIMER_HEADER = {
    "native": "timer_native.h",
    "gettimeofday": "timer_gettimeofday.h",
    "omp_get_wtime": "timer_omp_get_wtime.h",
    "mpi_wtime": "timer_mpi_wtime.h",
    "clock_gettime": "timer_clock_gettime.h",
    "papi_get_real_nsec": "timer_papi_real_nsec.h",
    "rdtsc": "timer_rdtsc.h",
    "rdtscp": "timer_rdtscp.h",
    "rdtscp_lfence": "timer_rdtscp_lfence.h",
    "tsc_asym": "timer_tsc_asym.h",
    "cntvct_el0": "timer_cntvct.h",
    "cntvct_el0_dmb": "timer_cntvct_dmb.h",
}

TIMER_PREFIX = {
    "native": "tacvar_timer_native",
    "gettimeofday": "tacvar_timer_gettimeofday",
    "omp_get_wtime": "tacvar_timer_omp_get_wtime",
    "mpi_wtime": "tacvar_timer_mpi_wtime",
    "clock_gettime": "tacvar_timer_clock_gettime",
    "papi_get_real_nsec": "tacvar_timer_papi_real_nsec",
    "rdtsc": "tacvar_timer_rdtsc",
    "rdtscp": "tacvar_timer_rdtscp",
    "rdtscp_lfence": "tacvar_timer_rdtscp_lfence",
    "tsc_asym": "tacvar_timer_tsc_asym",
    "cntvct_el0": "tacvar_timer_cntvct",
    "cntvct_el0_dmb": "tacvar_timer_cntvct_dmb",
}


def parse_conf(path: Path) -> dict[str, str]:
    cfg: dict[str, str] = {}
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if "=" not in s:
            raise SystemExit(f"bad conf line: {line}")
        k, v = s.split("=", 1)
        cfg[k.strip()] = v.strip().strip('"').strip("'")
    return cfg


def detect_arch(compiler: str) -> str:
    try:
        out = subprocess.check_output(
            [compiler, "-dM", "-E", "-"],
            input=b"",
            stderr=subprocess.DEVNULL,
        ).decode()
    except Exception as e:
        raise SystemExit(f"arch detect failed with {compiler}: {e}")
    if "__x86_64__" in out or "__amd64__" in out:
        return "x86_64"
    if "__aarch64__" in out:
        return "aarch64"
    raise SystemExit(f"unsupported target arch from {compiler}")


def parse_tf_mode(s: str) -> bool:
    v = (s or "OFF").strip().upper()
    if v in ("ON", "1", "TRUE", "YES"):
        return True
    if v in ("OFF", "0", "FALSE", "NO", ""):
        return False
    raise SystemExit(f"unknown TACVAR_TF_SAMPLING_MODE={s}")


def resolve_tf_root(raw: str, conf_dir: Path) -> str:
    if not raw:
        return ""
    p = Path(raw)
    if not p.is_absolute():
        p = (conf_dir / raw).resolve()
    else:
        p = p.resolve()
    return str(p)


def c_string(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def split_names(s: str) -> list[str]:
    if not s:
        return []
    return [p.strip() for p in s.split(",") if p.strip()]


def resolve_native(timer: str, consumer: str) -> str:
    if timer != "native":
        return timer
    if consumer == "npb-mpi":
        return "mpi_wtime"
    if consumer == "npb-omp":
        return "omp_get_wtime"
    if consumer == "lmbench":
        return "gettimeofday"
    raise SystemExit(f"native timer unsupported for consumer {consumer}")


def validate(cfg: dict[str, str], consumer: str, arch: str) -> tuple[str, str, list[str], int, int]:
    timer = cfg.get("TACVAR_TIMER", "native")
    counter = cfg.get("TACVAR_COUNTER_BACKEND", "none")
    names = split_names(cfg.get("TACVAR_COUNTER_NAMES", ""))
    nstp = int(cfg.get("TACVAR_NSTP", "0") or "0")
    count = int(cfg.get("TACVAR_COUNTER_COUNT", str(len(names))) or "0")

    if timer not in TIMERS_ALL:
        raise SystemExit(f"unknown TACVAR_TIMER={timer}")
    if counter not in COUNTERS:
        raise SystemExit(f"unknown TACVAR_COUNTER_BACKEND={counter}")
    if consumer not in CONSUMERS:
        raise SystemExit(f"unknown consumer={consumer}")

    # Consumer-restricted timers
    if timer == "mpi_wtime" and consumer != "npb-mpi":
        raise SystemExit("mpi_wtime only allowed for npb-mpi")
    if timer == "omp_get_wtime" and consumer != "npb-omp":
        raise SystemExit("omp_get_wtime only allowed for npb-omp")
    if timer == "gettimeofday" and consumer not in ("lmbench", "npb-omp"):
        # allow gettimeofday on lmbench; NPB-OMP may use via native fallback only
        if consumer != "lmbench":
            raise SystemExit("gettimeofday only allowed for lmbench (or via native on npb-omp)")

    resolved = resolve_native(timer, consumer)
    if resolved in TIMERS_X86 and arch != "x86_64":
        raise SystemExit(f"{resolved} requires x86_64")
    if resolved in TIMERS_ARM and arch != "aarch64":
        raise SystemExit(f"{resolved} requires aarch64")
    if resolved in TIMERS_ARM and nstp <= 0:
        raise SystemExit(f"{resolved} requires TACVAR_NSTP > 0")

    if count != len(names):
        raise SystemExit(
            f"TACVAR_COUNTER_COUNT={count} != len(TACVAR_COUNTER_NAMES)={len(names)}"
        )
    if count < 0 or count > 12:
        raise SystemExit("TACVAR_COUNTER_COUNT must be 0..12")
    if counter == "none" and count != 0:
        raise SystemExit("COUNTER_BACKEND=none requires COUNTER_COUNT=0")
    if counter != "none" and count == 0:
        raise SystemExit(f"{counter} requires COUNTER_COUNT > 0")
    if counter == "asm" and arch not in ("x86_64", "aarch64"):
        raise SystemExit("asm counters require x86_64 or aarch64")

    return timer, counter, names, count, nstp


def emit_header(
    out_h: Path,
    consumer: str,
    arch: str,
    timer: str,
    counter: str,
    names: list[str],
    count: int,
    nstp: int,
    output_root: str,
    measure_root: Path,
    tf_on: bool = False,
    tf_root: str = "",
    tf_nspg: float = 0.0,
    tf_reg: int = 0,
    tf_loc: int = 0,
) -> None:
    resolved = resolve_native(timer, consumer)
    prefix = TIMER_PREFIX[timer if timer == "native" else resolved]
    hdr = TIMER_HEADER[timer if timer == "native" else resolved]

    if counter == "none":
        c_hdr = "counter_none.h"
        c_prefix = "tacvar_counter_none"
    elif counter == "perf_event_open":
        c_hdr = "counter_perf_event.h"
        c_prefix = "tacvar_counter_perf_event"
    elif counter == "papi_read":
        c_hdr = "counter_papi.h"
        c_prefix = "tacvar_counter_papi"
    elif counter == "asm" and arch == "x86_64":
        c_hdr = "counter_rdpmc_x86.h"
        c_prefix = "tacvar_counter_rdpmc_x86"
    elif counter == "asm" and arch == "aarch64":
        c_hdr = "counter_pmu_arm.h"
        c_prefix = "tacvar_counter_pmu_arm"
    else:
        raise SystemExit(f"cannot map counter {counter}/{arch}")

    names_c = ", ".join(f'"{n}"' for n in names) if names else ""
    display_timer = timer if timer != "native" else f"native({resolved})"

    lines = [
        "/* Auto-generated by gen_config.py — do not edit. */",
        "#ifndef TACVAR_GENERATED_CONFIG_H_INCLUDED",
        "#define TACVAR_GENERATED_CONFIG_H_INCLUDED",
        "",
        f'#define TACVAR_CONSUMER_NPB_MPI 1',
        f'#define TACVAR_CONSUMER_NPB_OMP 2',
        f'#define TACVAR_CONSUMER_LMBENCH 3',
        f"#define TACVAR_CONSUMER {CONSUMERS[consumer]}",
        f'#define TACVAR_ARCH_STR "{arch}"',
        f'#define TACVAR_TIMER_NAME "{display_timer}"',
        f'#define TACVAR_TIMER_RESOLVED "{resolved}"',
        f'#define TACVAR_COUNTER_BACKEND_NAME "{counter}"',
        f"#define TACVAR_COUNTER_COUNT {count}",
        f"#define TACVAR_NSTP {nstp}",
        f'#define TACVAR_OUTPUT_ROOT_DEFAULT "{output_root}"',
        f"#define TACVAR_TF_SAMPLING {1 if tf_on else 0}",
        f'#define TACVAR_TF_DATA_ROOT "{c_string(tf_root)}"',
        f"#define TACVAR_TF_NSPG {tf_nspg}",
        f"#define TACVAR_TF_REG_ID {tf_reg}",
        f"#define TACVAR_TF_LOC_ID {tf_loc}",
        f"#define TACVAR_TF_TIMER_IS_TICK {1 if resolved in TIMERS_TICK else 0}",
        f"#define TACVAR_TF_REG_{tf_reg}_LOC_{tf_loc} 1",
        "#define TACVAR_TF_SIDE_OFF 0",
        "#define TACVAR_TF_SIDE_TFS 1",
        "#define TACVAR_TF_SIDE_TFE 2",
        "#ifndef TACVAR_TF_SIDE",
        "#define TACVAR_TF_SIDE TACVAR_TF_SIDE_OFF",
        "#endif",
        "",
    ]
    if names:
        lines.append(
            f"static const char *const TACVAR_COUNTER_NAME_LIST[{count}] = {{ {names_c} }};"
        )
    else:
        lines.append("/* no counter names */")
        lines.append(
            "static const char *const *const TACVAR_COUNTER_NAME_LIST = (const char *const *)0;"
        )

    lines += [
        "",
        f'#include "{hdr}"',
        f'#include "{c_hdr}"',
        "",
        f"#define TACVAR_TIMER_INIT()            {prefix}_init()",
        f"#define TACVAR_TIMER_FINI()            {prefix}_fini()",
        f"#define TACVAR_TIMER_BEGIN(raw)        {prefix}_begin(&(raw))",
        f"#define TACVAR_TIMER_END(raw)          {prefix}_end(&(raw))",
        f"#define TACVAR_TIMER_DELTA_NS(b, e)    {prefix}_delta_ns((b), (e))",
        "",
        f"#define TACVAR_COUNTER_INIT()          {c_prefix}_init(TACVAR_COUNTER_NAME_LIST, TACVAR_COUNTER_COUNT)",
        f"#define TACVAR_COUNTER_FINI()          {c_prefix}_fini()",
        f"#define TACVAR_COUNTER_READ(values)    {c_prefix}_read(values)",
        "",
        "#endif /* TACVAR_GENERATED_CONFIG_H_INCLUDED */",
        "",
    ]
    out_h.write_text("\n".join(lines))


def emit_mk(
    out_mk: Path,
    consumer: str,
    arch: str,
    timer: str,
    counter: str,
    papi_inc: str,
    papi_lib: str,
    measure_root: Path,
    tf_on: bool = False,
    tf_root: str = "",
    tf_nspg: float = 0.0,
    tf_reg: int = 0,
    tf_loc: int = 0,
) -> None:
    resolved = resolve_native(timer, consumer)
    objs = ["tacvar_measure.o", "tacvar_csv.o", "tacvar_tf.o"]
    cflags = [
        f"-I{measure_root}",
        f"-I{measure_root / 'include'}",
        f"-I{measure_root / 'timers'}",
        f"-I{measure_root / 'counters'}",
        f"-I{measure_root / 'events'}",
        f"-I{measure_root / 'gauges'}",
        "-DTACVAR_HAS_GENERATED_CONFIG=1",
    ]
    ldflags: list[str] = []
    need_tsc = resolved in TIMERS_X86 or (tf_on and arch == "x86_64")
    if need_tsc and "timer_util.o" not in objs:
        objs.append("timer_util.o")
    if counter == "perf_event_open":
        objs.append("counter_perf_event.o")
    elif counter == "papi_read":
        objs.append("counter_papi.o")
        if papi_inc:
            cflags.append(f"-I{papi_inc}")
        if papi_lib:
            ldflags += [f"-L{papi_lib}", "-lpapi", f"-Wl,-rpath,{papi_lib}"]
        else:
            ldflags.append("-lpapi")
    elif counter == "asm" and arch == "x86_64":
        objs += ["counter_rdpmc_x86.o", "x86_events.o"]
    elif counter == "asm" and arch == "aarch64":
        objs += ["counter_pmu_arm.o", "armv8_events.o"]

    if resolved == "papi_get_real_nsec" or timer == "papi_get_real_nsec":
        if papi_inc:
            cflags.append(f"-I{papi_inc}")
        if "-lpapi" not in " ".join(ldflags):
            if papi_lib:
                ldflags += [f"-L{papi_lib}", "-lpapi", f"-Wl,-rpath,{papi_lib}"]
            else:
                ldflags.append("-lpapi")

    if resolved == "mpi_wtime" or (timer == "native" and consumer == "npb-mpi"):
        pass  # MPI via mpicc

    lines = [
        "# Auto-generated by gen_config.py — do not edit.",
        f"TACVAR_ARCH := {arch}",
        f"TACVAR_CONSUMER := {consumer}",
        f"TACVAR_TIMER := {timer}",
        f"TACVAR_TIMER_RESOLVED := {resolved}",
        f"TACVAR_COUNTER_BACKEND := {counter}",
        f"TACVAR_MEASURE_OBJS := {' '.join(objs)}",
        f"TACVAR_MEASURE_CFLAGS := {' '.join(cflags)}",
        f"TACVAR_MEASURE_LDFLAGS := {' '.join(ldflags)}",
        f"TACVAR_MEASURE_ROOT := {measure_root}",
        f"TACVAR_TF_SAMPLING_MODE := {'ON' if tf_on else 'OFF'}",
        f"TACVAR_TF_DATA_ROOT := {tf_root}",
        f"TACVAR_TF_NSPG := {tf_nspg}",
        f"TACVAR_TF_REG_ID := {tf_reg}",
        f"TACVAR_TF_LOC_ID := {tf_loc}",
        "",
    ]
    out_mk.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--conf", required=True)
    ap.add_argument("--compiler", required=True)
    ap.add_argument("--consumer", required=True, choices=list(CONSUMERS))
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--measure-root", required=True)
    ap.add_argument("--papi-inc", default="")
    ap.add_argument("--papi-lib", default="")
    args = ap.parse_args()

    conf = Path(args.conf)
    outdir = Path(args.outdir)
    measure_root = Path(args.measure_root).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cfg = parse_conf(conf)
    arch = detect_arch(args.compiler)
    timer, counter, names, count, nstp = validate(cfg, args.consumer, arch)
    output_root = cfg.get("TACVAR_OUTPUT_ROOT", ".")
    tf_on = parse_tf_mode(cfg.get("TACVAR_TF_SAMPLING_MODE", "OFF"))
    tf_root = resolve_tf_root(cfg.get("TACVAR_TF_DATA_ROOT", ""), conf.parent)
    try:
        tf_nspg = float(cfg.get("TACVAR_TF_NSPG", "0") or "0")
    except ValueError:
        raise SystemExit("TACVAR_TF_NSPG must be a number")
    try:
        tf_reg = int(cfg.get("TACVAR_TF_REG_ID", "0") or "0")
        tf_loc = int(cfg.get("TACVAR_TF_LOC_ID", "0") or "0")
    except ValueError:
        raise SystemExit("TACVAR_TF_REG_ID and TACVAR_TF_LOC_ID must be integers")
    if tf_reg < 0 or tf_loc < 0:
        raise SystemExit("TACVAR_TF_REG_ID and TACVAR_TF_LOC_ID must be >= 0")
    if tf_on and arch == "aarch64" and nstp <= 0:
        raise SystemExit(
            "TACVAR_TF_SAMPLING_MODE=ON on aarch64 requires TACVAR_NSTP > 0 "
            "(tick-to-ns for unfenced cntvct)"
        )

    # PAPI path required when selected
    if counter == "papi_read" or timer == "papi_get_real_nsec":
        if not args.papi_inc and not os.environ.get("PAPI_INC"):
            # still allow system headers
            pass

    emit_header(
        outdir / "tacvar_generated_config.h",
        args.consumer,
        arch,
        timer,
        counter,
        names,
        count,
        nstp,
        output_root,
        measure_root,
        tf_on,
        tf_root,
        tf_nspg,
        tf_reg,
        tf_loc,
    )
    emit_mk(
        outdir / "tacvar_generated_config.mk",
        args.consumer,
        arch,
        timer,
        counter,
        args.papi_inc or os.environ.get("PAPI_INC", ""),
        args.papi_lib or os.environ.get("PAPI_LIB", ""),
        measure_root,
        tf_on,
        tf_root,
        tf_nspg,
        tf_reg,
        tf_loc,
    )
    extra = f" tf={('ON' if tf_on else 'OFF')}"
    print(
        f"generated config for {args.consumer}/{arch} timer={timer} "
        f"counter={counter}{extra}"
    )


if __name__ == "__main__":
    main()
