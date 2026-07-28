#!/usr/bin/env bash
# Pairwise NPB-MPI timer/counter matrix (MG Class S). Not Cartesian.
# Classifies PASS / SKIP (host limitation) / FAIL.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$ROOT/../.." && pwd)"
MPI_DIR="$ROOT/NPB3.4-MPI"
ARCH=$(uname -m)
LOG="${TACVAR_MATRIX_LOG:-/tmp/npb_backend_matrix_${ARCH}.log}"
: > "$LOG"

export PATH="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.7}/bin:${PATH}"
export LD_LIBRARY_PATH="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.7}/lib:${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}/lib:${LD_LIBRARY_PATH:-}"
NSTP=0
if [[ "$ARCH" == "aarch64" ]]; then
  export PATH="/home/hpckey/01-App/openmpi-5.0.8/bin:${PATH}"
  export LD_LIBRARY_PATH="/home/hpckey/01-App/openmpi-5.0.8/lib:/home/hpckey/01-App/papi/lib:${LD_LIBRARY_PATH:-}"
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi}"
  export MPI_HOME="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.8}"
  # Derive ns/tick from CNTFRQ_EL0 when possible.
  if [[ -r /sys/devices/system/cpu/cpu0/cpu_capacity ]] || true; then
    frq=$(awk '/CNTFRQ|cntfrq/ {print $3; exit}' /proc/cpuinfo 2>/dev/null || true)
    if [[ -z "${frq:-}" ]]; then
      # Common HiSilicon CNTVCT frequency is 100 MHz → 10 ns/tick
      frq=100000000
    fi
    if [[ "$frq" =~ ^[0-9]+$ ]] && [[ "$frq" -gt 0 ]]; then
      NSTP=$(( (1000000000 + frq/2) / frq ))
    else
      NSTP=10
    fi
  fi
  echo "ARM TACVAR_NSTP=$NSTP (from CNTVCT freq assumption ${frq:-unknown})" | tee -a "$LOG"
else
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}"
  export MPI_HOME="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.7}"
fi

pass=0
skip=0
fail=0

apply_conf() {
  local timer="$1" counter="$2" names="$3" count="$4" nstp="$5"
  cat > "$MPI_DIR/tacvar.conf.matrix.bak" <<EOF
# backup marker
EOF
  cp "$MPI_DIR/tacvar.conf" "$MPI_DIR/tacvar.conf.matrix.bak"
  cat > "$MPI_DIR/tacvar.conf" <<EOF
TACVAR_TIMER=$timer
TACVAR_NSTP=$nstp
TACVAR_COUNTER_BACKEND=$counter
TACVAR_COUNTER_COUNT=$count
TACVAR_COUNTER_NAMES=$names
TACVAR_OUTPUT_ROOT=.
TACVAR_ENABLE_PER_STEP_TIMING=1
EOF
}

restore_conf() {
  if [[ -f "$MPI_DIR/tacvar.conf.matrix.bak" ]]; then
    mv -f "$MPI_DIR/tacvar.conf.matrix.bak" "$MPI_DIR/tacvar.conf"
    touch "$MPI_DIR/tacvar.conf"
  fi
}

force_rebuild() {
  local arch
  arch=$(uname -m)
  rm -rf "$MPI_DIR/config/tacvar_build" "$MPI_DIR/config/tacvar_build_${arch}"
  rm -f "$MPI_DIR/common/"*.o "$MPI_DIR/sys/setparams" "$MPI_DIR/bin/"*.x
  find "$MPI_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.o' -delete 2>/dev/null || true
}

is_skip_output() {
  local out="$1"
  # Host limitations only — not generic compile missing-headers from our Makefiles.
  grep -qiE 'Permission denied|Operation not permitted|not supported|PAPI_.*disabled|Component .* disabled|cannot open ph_enable_pmu|Event does not exist|Event.*not available|EPERM|ENOTSUP|perf_event_open failed|Failed to open|No such device|kmod' "$out" 2>/dev/null
}

check_run() {
  local label="$1" out="$2" expect_counters="$3"
  local dir csv steps totals nz

  if ! grep -q 'SUCCESSFUL' "$out"; then
    if is_skip_output "$out"; then
      echo "SKIP: $label (host limitation)" | tee -a "$LOG"
      skip=$((skip+1))
      return 0
    fi
    echo "FAIL: $label (verification)" | tee -a "$LOG"
    fail=$((fail+1))
    return 1
  fi

  dir=$(ls -dt "$MPI_DIR"/data_????????T?????? 2>/dev/null | head -1 || true)
  if [[ -z "$dir" || ! -f "$dir/region_info.csv" ]]; then
    echo "FAIL: $label (missing data/region_info)" | tee -a "$LOG"
    fail=$((fail+1))
    return 1
  fi

  for csv in "$dir"/npb-mpi_mg_*.csv; do
    [[ -f "$csv" ]] || { echo "FAIL: $label (no event csv)"; fail=$((fail+1)); return 1; }
    steps=$(awk -F, 'NR>1 && $6+0==1000 {n++} END{print n+0}' "$csv")
    totals=$(awk -F, 'NR>1 && $6+0!=1000 {n++} END{print n+0}' "$csv")
    if [[ "$steps" != "4" || "$totals" != "1" ]]; then
      echo "FAIL: $label (events steps=$steps totals=$totals)" | tee -a "$LOG"
      fail=$((fail+1))
      return 1
    fi
    if [[ "$expect_counters" == "1" ]]; then
      nz=$(awk -F, '
        NR==1 { for(i=1;i<=NF;i++) if($i ~ /_delta$/) d[++nd]=i; next }
        { for(j=1;j<=nd;j++) if($(d[j])+0 != 0) nz=1 }
        END { print (nd>=2 && nz) ? 1 : 0 }
      ' "$csv")
      if [[ "$nz" != "1" ]]; then
        if is_skip_output "$out"; then
          echo "SKIP: $label (zero counters / host)" | tee -a "$LOG"
          skip=$((skip+1))
          return 0
        fi
        echo "FAIL: $label (zero/missing counter deltas)" | tee -a "$LOG"
        fail=$((fail+1))
        return 1
      fi
    fi
  done
  echo "PASS: $label" | tee -a "$LOG"
  pass=$((pass+1))
  return 0
}

run_combo() {
  local timer="$1" counter="$2" names="$3" count="$4" nstp="$5"
  local label="${timer}+${counter}"
  local out="/tmp/npb_matrix_${ARCH}_${timer}_${counter}.out"
  local expect_c=0
  [[ "$count" -gt 0 ]] && expect_c=1

  echo "=== $label ===" | tee -a "$LOG"
  apply_conf "$timer" "$counter" "$names" "$count" "$nstp"
  force_rebuild
  set +e
  (cd "$MPI_DIR" && make MG CLASS=S) >"$out" 2>&1
  rc=$?
  if [[ $rc -ne 0 ]]; then
    if is_skip_output "$out"; then
      echo "SKIP: $label (build host limitation)" | tee -a "$LOG"
      skip=$((skip+1))
      set -e
      restore_conf
      return 0
    fi
    echo "FAIL: $label (compile)" | tee -a "$LOG"
    tail -40 "$out" | tee -a "$LOG"
    fail=$((fail+1))
    set -e
    restore_conf
    return 0
  fi
  rm -rf "$MPI_DIR"/data_????????T??????
  (cd "$MPI_DIR" && mpirun -np 4 --bind-to core ./bin/mg.S.x) >>"$out" 2>&1
  rc=$?
  set -e
  if [[ $rc -ne 0 ]] && is_skip_output "$out"; then
    echo "SKIP: $label (runtime host limitation)" | tee -a "$LOG"
    skip=$((skip+1))
    restore_conf
    return 0
  fi
  check_run "$label" "$out" "$expect_c" || true
  restore_conf
}

trap 'restore_conf; force_rebuild' EXIT

echo "NPB-MPI backend matrix arch=$ARCH log=$LOG" | tee -a "$LOG"

if [[ "$ARCH" == "x86_64" ]]; then
  # Pairwise: each timer and each reader at least once
  run_combo native none "" 0 0
  run_combo mpi_wtime perf_event_open "cpu-cycles,instructions" 2 0
  run_combo clock_gettime papi_read "PAPI_TOT_CYC,PAPI_TOT_INS" 2 0
  run_combo papi_get_real_nsec none "" 0 0
  run_combo rdtsc none "" 0 0
  run_combo rdtscp none "" 0 0
  run_combo rdtscp_lfence none "" 0 0
  run_combo tsc_asym asm "INST_RETIRED.ANY_P,CPU_CLK_UNHALTED.THREAD" 2 0
elif [[ "$ARCH" == "aarch64" ]]; then
  run_combo native none "" 0 0
  run_combo mpi_wtime perf_event_open "cpu-cycles,instructions" 2 0
  run_combo clock_gettime papi_read "PAPI_TOT_CYC,PAPI_TOT_INS" 2 0
  run_combo papi_get_real_nsec none "" 0 0
  run_combo cntvct_el0 none "" 0 "$NSTP"
  run_combo cntvct_el0_dmb asm "CPU_CYCLES,INST_RETIRED" 2 "$NSTP"
else
  echo "Unsupported arch $ARCH" | tee -a "$LOG"
  exit 2
fi

echo "matrix summary: pass=$pass skip=$skip fail=$fail" | tee -a "$LOG"
[[ "$fail" -eq 0 ]]
