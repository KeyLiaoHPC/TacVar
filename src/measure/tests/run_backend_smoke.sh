#!/usr/bin/env bash
# Backend smoke: build/read each available timer/counter on this arch.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ARCH=$(uname -m)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}"

build_one() {
  local timer="$1" counter="$2" names="$3" count="$4" nstp="$5"
  local dir="$TMP/${timer}_${counter}"
  mkdir -p "$dir"
  cat > "$dir/tacvar.conf" <<EOF
TACVAR_TIMER=$timer
TACVAR_NSTP=$nstp
TACVAR_COUNTER_BACKEND=$counter
TACVAR_COUNTER_COUNT=$count
TACVAR_COUNTER_NAMES=$names
TACVAR_OUTPUT_ROOT=$dir
EOF
  make -C "$ROOT" CONF="$dir/tacvar.conf" CONSUMER=lmbench OUTDIR="$dir/build" \
    CC=gcc PAPI_INC="$PAPI_HOME/include" PAPI_LIB="$PAPI_HOME/lib" >/dev/null
  echo "built $timer + $counter"
}

echo "arch=$ARCH"
build_one native none "" 0 0
build_one clock_gettime none "" 0 0
build_one gettimeofday none "" 0 0
if [[ -f "$PAPI_HOME/include/papi.h" ]]; then
  build_one papi_get_real_nsec none "" 0 0 || echo "skip papi timer"
fi
if [[ "$ARCH" == "x86_64" ]]; then
  build_one rdtsc none "" 0 0
  build_one rdtscp none "" 0 0
  build_one rdtscp_lfence none "" 0 0
  build_one tsc_asym none "" 0 0
  build_one clock_gettime perf_event_open "cpu-cycles,instructions" 2 0 || echo "skip perf"
  if [[ -d /sys/module/ph_enable_pmu ]]; then
    build_one tsc_asym asm "INST_RETIRED.ANY_P,CPU_CLK_UNHALTED.THREAD" 2 0 || echo "skip rdpmc"
  else
    echo "skip asm/rdpmc (kmod not loaded)"
  fi
elif [[ "$ARCH" == "aarch64" ]]; then
  build_one cntvct_el0 none "" 0 10
  build_one cntvct_el0_dmb none "" 0 10
  build_one clock_gettime perf_event_open "cpu-cycles,instructions" 2 0 || echo "skip perf"
  if [[ -d /sys/module/ph_enable_pmu ]]; then
    build_one cntvct_el0 asm "CPU_CYCLES,INST_RETIRED" 2 10 || echo "skip arm asm"
  else
    echo "skip arm asm (kmod not loaded)"
  fi
fi
echo "backend smoke done"
