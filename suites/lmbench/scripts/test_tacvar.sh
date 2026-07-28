#!/usr/bin/env bash
# lmbench TacVar build + smoke (single-core lat_syscall null).
set -euo pipefail
LMBENCH="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$LMBENCH/../.." && pwd)"
MODE="${1:---build-only}"
ARCH=$(uname -m)
PROTECT="$REPO/src/measure/tests/protected_sources.txt"
SNAPSHOT=$(mktemp)
trap 'rm -f "$SNAPSHOT"' EXIT

if [[ "$ARCH" == "aarch64" ]]; then
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.2.0b2}"
else
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}"
fi
export LD_LIBRARY_PATH="${PAPI_HOME}/lib:${LD_LIBRARY_PATH:-}"

# scripts/os must be invoked from src/ so it finds scripts/gnu-os
OS=$(cd "$LMBENCH/src" && ../scripts/os)
BINDIR="$LMBENCH/bin/$OS"
RELBIN="../bin/${OS}"

snapshot_protected() {
  (cd "$REPO" && while IFS= read -r p; do
    [[ -z "$p" || "$p" =~ ^# ]] && continue
    [[ -f "$p" ]] && sha256sum "$p"
  done < "$PROTECT") > "$SNAPSHOT"
}

check_protected() {
  local now
  now=$(mktemp)
  (cd "$REPO" && while IFS= read -r p; do
    [[ -z "$p" || "$p" =~ ^# ]] && continue
    [[ -f "$p" ]] && sha256sum "$p"
  done < "$PROTECT") > "$now"
  if ! diff -q "$SNAPSHOT" "$now" >/dev/null; then
    echo "ERROR: protected lmbench sources were modified"
    diff -u "$SNAPSHOT" "$now" || true
    rm -f "$now"
    exit 1
  fi
  rm -f "$now"
  echo "protected sources unchanged"
}

check_csv() {
  local cwd="$1"
  local dir
  dir=$(ls -dt "$cwd"/data_????????T?????? 2>/dev/null | head -1 || true)
  [[ -n "$dir" ]] || { echo "ERROR: no data_* under $cwd"; exit 1; }
  local csv
  csv=$(find "$dir" -name '*.csv' | head -1)
  [[ -n "$csv" ]] || { echo "ERROR: no CSV"; exit 1; }
  head -1 "$csv" | grep -q elapsed_ns
  awk -F, 'NR>1 && $10+0 < 0 {exit 1}' "$csv"
  echo "CSV OK: $csv"
}

apply_temp_conf() {
  local timer="$1" counter="$2" names="$3" count="$4" nstp="$5"
  local bak="$LMBENCH/tacvar.conf.smoke.bak"
  cp "$LMBENCH/tacvar.conf" "$bak"
  cat > "$LMBENCH/tacvar.conf" <<EOF
TACVAR_TIMER=$timer
TACVAR_NSTP=$nstp
TACVAR_COUNTER_BACKEND=$counter
TACVAR_COUNTER_COUNT=$count
TACVAR_COUNTER_NAMES=$names
TACVAR_OUTPUT_ROOT=.
EOF
  echo "$bak"
}

restore_conf() {
  local bak="$1"
  if [[ -f "$bak" ]]; then
    mv -f "$bak" "$LMBENCH/tacvar.conf"
  fi
  return 0
}

build_targets() {
  # Relative $O paths required (scripts/build sets O=../bin/$OS).
  # Force rebuild so conf changes take effect.
  rm -f "$BINDIR/lib_timing.o" "$BINDIR/tacvar_lmbench.o" "$BINDIR/lmbench.a" \
        "$BINDIR/lat_syscall" "$BINDIR/timing_o" "$BINDIR/enough"
  rm -rf "$BINDIR/tacvar_build"
  (cd "$LMBENCH/src" && ../scripts/build \
    "${RELBIN}/lmbench.a" "${RELBIN}/timing_o" "${RELBIN}/enough" "${RELBIN}/lat_syscall")
}

snapshot_protected

echo "=== lmbench build (native+none) ==="
build_targets
[[ -x "$BINDIR/lat_syscall" ]]
[[ -f "$BINDIR/lmbench.a" ]] || [[ -f "$LMBENCH/src/lmbench.a" ]] || true
if [[ -f "$BINDIR/lib_timing.o" ]] || [[ -f "$LMBENCH/src/lib_timing.o" ]]; then
  obj=${BINDIR}/lib_timing.o
  [[ -f "$obj" ]] || obj="$LMBENCH/src/lib_timing.o"
  if objdump -d "$obj" 2>/dev/null | grep -E 'call\s+\*%|blr\s+x' >/dev/null; then
    echo "NOTE: indirect call in $obj"
  else
    echo "hotpath OK: $obj"
  fi
fi

if [[ "$MODE" == "--build-only" ]]; then
  check_protected
  echo "lmbench build-only OK"
  exit 0
fi

if [[ "$MODE" != "--run-smoke" ]]; then
  echo "Usage: $0 [--build-only|--run-smoke]"
  exit 2
fi

echo "=== lmbench smoke: lat_syscall null ==="
rm -rf "$LMBENCH"/data_????????T??????
SMOKE_OUT=$(mktemp)
( cd "$LMBENCH" && taskset -c 0 "$BINDIR/lat_syscall" null ) | tee "$SMOKE_OUT"
if ! grep -Eiq 'Simple syscall|microseconds|usecs|null' "$SMOKE_OUT"; then
  echo "WARN: unexpected lat_syscall output format"
  cat "$SMOKE_OUT"
fi
rm -f "$SMOKE_OUT"
check_csv "$LMBENCH"

if [[ "$ARCH" == "x86_64" ]]; then
  if [[ -d /sys/module/ph_enable_pmu ]]; then
    echo "=== lmbench combo: tsc_asym + asm ==="
    bak=$(apply_temp_conf tsc_asym asm "INST_RETIRED.ANY_P,CPU_CLK_UNHALTED.THREAD" 2 0)
    trap 'restore_conf "'"$bak"'"; rm -f '"$SNAPSHOT" EXIT
    build_targets
    rm -rf "$LMBENCH"/data_????????T??????
    ( cd "$LMBENCH" && taskset -c 0 "$BINDIR/lat_syscall" null ) | tee /tmp/lmb_combo.out
    check_csv "$LMBENCH"
    restore_conf "$bak"
    build_targets
    trap 'rm -f '"$SNAPSHOT" EXIT
  else
    echo "skip tsc_asym+asm (kmod not loaded); trying clock_gettime+perf"
    bak=$(apply_temp_conf clock_gettime perf_event_open "cpu-cycles,instructions" 2 0)
    trap 'restore_conf "'"$bak"'"; rm -f '"$SNAPSHOT" EXIT
    build_targets
    rm -rf "$LMBENCH"/data_????????T??????
    ( cd "$LMBENCH" && taskset -c 0 "$BINDIR/lat_syscall" null )
    check_csv "$LMBENCH"
    restore_conf "$bak"
    build_targets
    trap 'rm -f '"$SNAPSHOT" EXIT
  fi
elif [[ "$ARCH" == "aarch64" ]]; then
  NSTP="${TACVAR_NSTP_ARM:-10}"
  echo "=== lmbench combo: cntvct_el0 + perf_event_open (nstp=$NSTP) ==="
  bak=$(apply_temp_conf cntvct_el0 perf_event_open "cpu-cycles,instructions" 2 "$NSTP")
  trap 'restore_conf "'"$bak"'"; rm -f '"$SNAPSHOT" EXIT
  build_targets
  rm -rf "$LMBENCH"/data_????????T??????
  ( cd "$LMBENCH" && taskset -c 0 "$BINDIR/lat_syscall" null )
  check_csv "$LMBENCH"
  restore_conf "$bak"
  build_targets
  trap 'rm -f '"$SNAPSHOT" EXIT
fi

check_protected
echo "lmbench smoke OK"
exit 0
