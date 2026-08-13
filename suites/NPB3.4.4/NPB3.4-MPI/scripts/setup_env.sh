#!/usr/bin/env bash
# Prepare a node for TacVar NPB-MPI timing: toolchain, frequency, PMU access.
# Source or run:  source scripts/setup_env.sh   OR   bash scripts/setup_env.sh
set -euo pipefail

ARCH=$(uname -m)
HOST=$(hostname -s 2>/dev/null || hostname)
APP="${HOME}/01-App"

if [[ "$ARCH" == "aarch64" ]]; then
  export MPI_HOME="${MPI_HOME:-${APP}/openmpi-5.0.8}"
  export PAPI_HOME="${PAPI_HOME:-${APP}/papi-7.2.0b2}"
  export TACVAR_NSTP_ARM="${TACVAR_NSTP_ARM:-10}"
else
  export MPI_HOME="${MPI_HOME:-${APP}/openmpi-5.0.7}"
  export PAPI_HOME="${PAPI_HOME:-${APP}/papi-7.1.0}"
fi

if [[ -d "${MPI_HOME}/bin" ]]; then
  export PATH="${MPI_HOME}/bin:${PATH}"
fi
_ld=""
[[ -d "${MPI_HOME}/lib" ]] && _ld="${MPI_HOME}/lib"
[[ -d "${PAPI_HOME}/lib" ]] && _ld="${_ld:+$_ld:}${PAPI_HOME}/lib"
if [[ -n "$_ld" ]]; then
  export LD_LIBRARY_PATH="${_ld}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

echo "setup_env: host=${HOST} arch=${ARCH}"
echo "setup_env: MPI_HOME=${MPI_HOME}"
echo "setup_env: PAPI_HOME=${PAPI_HOME}"
command -v mpirun >/dev/null && echo "setup_env: mpirun=$(command -v mpirun)" || echo "setup_env: mpirun not on PATH"

# Fixed frequency / no turbo (best-effort; may need root).
if [[ -w /sys/devices/system/cpu/intel_pstate/no_turbo ]]; then
  echo 1 >/sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || true
  echo "setup_env: intel_pstate no_turbo=$(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo n/a)"
fi
if [[ -d /sys/devices/system/cpu/cpu0/cpufreq ]]; then
  gov=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)
  echo "setup_env: cpufreq governor=${gov} (set to performance as root if needed)"
fi

if [[ "$ARCH" == "aarch64" ]]; then
  echo "setup_env: TACVAR_NSTP_ARM=${TACVAR_NSTP_ARM} (KunPeng CNTVCT is typically 10 ns/tick)"
  if [[ -r /sys/devices/system/cpu/cpu0/cpu_capacity ]]; then
    echo "setup_env: cpu0 capacity=$(cat /sys/devices/system/cpu/cpu0/cpu_capacity)"
  fi
fi

if [[ -d /sys/module/ph_enable_pmu ]]; then
  echo "setup_env: ph_enable_pmu loaded"
else
  echo "setup_env: ph_enable_pmu not loaded (needed for asm/RDPMC counters)"
fi

if [[ -r /proc/sys/kernel/perf_event_paranoid ]]; then
  echo "setup_env: perf_event_paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)"
fi
if [[ -r /sys/bus/event_source/devices/cpu/rdpmc ]]; then
  echo "setup_env: rdpmc=$(cat /sys/bus/event_source/devices/cpu/rdpmc)"
fi
