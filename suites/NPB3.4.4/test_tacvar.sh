#!/usr/bin/env bash
# NPB TacVar build + smoke (<=4 cores). Does not modify benchmark kernels.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$ROOT/../.." && pwd)"
MODE="${1:---build-only}"
ARCH=$(uname -m)
PROTECT="$REPO/src/measure/tests/protected_sources.txt"
SNAPSHOT=$(mktemp)
trap 'rm -f "$SNAPSHOT"' EXIT

export PATH="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.7}/bin:${PATH}"
export LD_LIBRARY_PATH="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.7}/lib:${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}/lib:${LD_LIBRARY_PATH:-}"
if [[ "$ARCH" == "aarch64" ]]; then
  export PATH="/home/hpckey/01-App/openmpi-5.0.8/bin:${PATH}"
  export LD_LIBRARY_PATH="/home/hpckey/01-App/openmpi-5.0.8/lib:/home/hpckey/01-App/papi/lib:${LD_LIBRARY_PATH:-}"
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi}"
  export MPI_HOME="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.8}"
else
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}"
  export MPI_HOME="${MPI_HOME:-/home/hpckey/01-App/openmpi-5.0.7}"
fi

snapshot_protected() {
  (cd "$REPO" && while IFS= read -r p; do
    [[ -z "$p" || "$p" =~ ^# ]] && continue
    if [[ -f "$p" ]]; then
      sha256sum "$p"
    elif [[ -d "$p" ]]; then
      find "$p" -type f \( -name '*.f' -o -name '*.f90' -o -name '*.c' -o -name '*.h' \) -print0 \
        | sort -z | xargs -0 sha256sum 2>/dev/null || true
    fi
  done < "$PROTECT") > "$SNAPSHOT"
}

check_protected() {
  local now
  now=$(mktemp)
  snapshot_protected
  # re-snapshot into $now
  (cd "$REPO" && while IFS= read -r p; do
    [[ -z "$p" || "$p" =~ ^# ]] && continue
    if [[ -f "$p" ]]; then sha256sum "$p"
    elif [[ -d "$p" ]]; then
      find "$p" -type f \( -name '*.f' -o -name '*.f90' -o -name '*.c' -o -name '*.h' \) -print0 \
        | sort -z | xargs -0 sha256sum 2>/dev/null || true
    fi
  done < "$PROTECT") > "$now"
  if ! diff -q "$SNAPSHOT" "$now" >/dev/null; then
    echo "ERROR: protected NPB sources were modified by the test"
    diff -u "$SNAPSHOT" "$now" || true
    rm -f "$now"
    exit 1
  fi
  rm -f "$now"
  echo "protected sources unchanged"
}

check_csv_dir() {
  local cwd="$1" expect_ranks="${2:-}"
  local dir
  dir=$(ls -dt "$cwd"/data_????????T?????? 2>/dev/null | head -1 || true)
  [[ -n "$dir" ]] || { echo "ERROR: no data_* under $cwd"; exit 1; }
  local csvs
  mapfile -t csvs < <(find "$dir" -name '*.csv' ! -name 'region_info.csv' | sort)
  [[ ${#csvs[@]} -ge 1 ]] || { echo "ERROR: no event CSV in $dir"; exit 1; }
  local hdr
  hdr=$(head -1 "${csvs[0]}")
  echo "$hdr" | grep -q 'elapsed_ns' || { echo "ERROR: bad CSV header"; exit 1; }
  echo "$hdr" | grep -q 'timer' || { echo "ERROR: missing timer column"; exit 1; }
  # non-negative elapsed
  awk -F, 'NR>1 && $10+0 < 0 { bad=1 } END { exit bad+0 }' "${csvs[@]}" \
    || { echo "ERROR: negative elapsed_ns"; exit 1; }
  if [[ -n "$expect_ranks" ]]; then
    local ranks
    ranks=$(awk -F, 'NR>1 {print $11}' "${csvs[@]}" | sort -nu | tr '\n' ',')
    echo "CSV ranks: $ranks (expect 0..$((expect_ranks-1)))"
  fi
  echo "CSV OK: $dir (${#csvs[@]} files)"
}

check_per_step_mg_cg() {
  local cwd="$1"
  local dir steps region_info
  dir=$(ls -dt "$cwd"/data_????????T?????? 2>/dev/null | head -1 || true)
  [[ -n "$dir" ]] || { echo "ERROR: no data_* under $cwd"; exit 1; }
  region_info="$dir/region_info.csv"
  [[ -f "$region_info" ]] || { echo "ERROR: missing region_info.csv"; exit 1; }
  grep -q ',1000,' "$region_info" || { echo "ERROR: region_info missing step region 1000"; exit 1; }
  # Each rank event file must have exactly expect_steps rows with region_id=1000
  local csv
  for csv in "$dir"/npb-mpi_*.csv; do
    [[ -f "$csv" ]] || continue
    [[ "$(basename "$csv")" == region_info.csv ]] && continue
    steps=$(awk -F, 'NR>1 && $6+0==1000 {n++} END{print n+0}' "$csv")
    local bench
    bench=$(awk -F, 'NR==2 {print $3; exit}' "$csv")
    local expect=0
    case "$bench" in
      mg) expect=4 ;;
      cg) expect=15 ;;
      *) echo "ERROR: unexpected bench $bench in $csv"; exit 1 ;;
    esac
    [[ "$steps" == "$expect" ]] || {
      echo "ERROR: $csv has $steps step rows, expect $expect"
      exit 1
    }
    # require two counter delta columns present and some non-zero
    local hdr
    hdr=$(head -1 "$csv")
    echo "$hdr" | grep -q '_delta' || { echo "ERROR: missing counter deltas in $csv"; exit 1; }
    awk -F, -v nsteps="$expect" '
      NR==1 {
        for(i=1;i<=NF;i++) if($i ~ /_delta$/) d[++nd]=i
        next
      }
      $6+0==1000 {
        for(j=1;j<=nd;j++) if($(d[j])+0 != 0) nz=1
      }
      END { exit (nd>=2 && nz) ? 0 : 1 }
    ' "$csv" || { echo "ERROR: need >=2 counter deltas with non-zero step values in $csv"; exit 1; }
  done
  echo "per-step OK: $dir"
}

apply_temp_conf() {
  local suite_dir="$1" timer="$2" counter="$3" names="$4" count="$5" nstp="$6"
  local per_step="${7:-0}"
  local bak="$suite_dir/tacvar.conf.smoke.bak"
  cp "$suite_dir/tacvar.conf" "$bak"
  cat > "$suite_dir/tacvar.conf" <<EOF
TACVAR_TIMER=$timer
TACVAR_NSTP=$nstp
TACVAR_COUNTER_BACKEND=$counter
TACVAR_COUNTER_COUNT=$count
TACVAR_COUNTER_NAMES=$names
TACVAR_OUTPUT_ROOT=.
TACVAR_ENABLE_PER_STEP_TIMING=$per_step
EOF
  echo "$bak"
}

restore_conf() {
  local bak="$1"
  local conf="${bak%.smoke.bak}"
  if [[ -f "$bak" ]]; then
    mv -f "$bak" "$conf"
    # Ensure make sees conf as newer than generated artifacts.
    touch "$conf"
  fi
}

force_rebuild_measure() {
  local suite_dir="$1"
  local arch
  arch=$(uname -m)
  rm -rf "$suite_dir/config/tacvar_build" "$suite_dir/config/tacvar_build_${arch}"
  rm -f "$suite_dir/common/"*.o
  # Host-arch binaries/objects (NFS shared with the other ISA)
  rm -f "$suite_dir/sys/setparams" "$suite_dir/bin/"*.x
  # Benchmark .o from the other arch (EM_AARCH64 vs x86_64)
  find "$suite_dir" -mindepth 2 -maxdepth 2 -type f -name '*.o' -delete 2>/dev/null || true
}

build_omp() {
  force_rebuild_measure "$ROOT/NPB3.4-OMP"
  (cd "$ROOT/NPB3.4-OMP" && make CG CLASS=S)
}

build_mpi() {
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  (cd "$ROOT/NPB3.4-MPI" && make IS CLASS=S)
}

hotpath_check() {
  local obj="$1"
  [[ -f "$obj" ]] || return 0
  if objdump -d "$obj" 2>/dev/null | grep -E 'call\s+\*%|jmp\s+\*%|blr\s+x' >/dev/null; then
    # Heuristic only: warn, do not fail (libraries may have unrelated indirect calls)
    echo "NOTE: indirect call/jump in $obj (review if timer hot path)"
  else
    echo "hotpath OK (no obvious indirect dispatch): $obj"
  fi
}

snapshot_protected

echo "=== NPB build (native+none from suite tacvar.conf) ==="
build_omp
build_mpi
hotpath_check "$ROOT/NPB3.4-OMP/common/c_timers.o"
hotpath_check "$ROOT/NPB3.4-MPI/common/c_timers.o"

if [[ "$MODE" == "--build-only" ]]; then
  check_protected
  echo "NPB build-only OK"
  exit 0
fi

if [[ "$MODE" != "--run-smoke" ]]; then
  echo "Usage: $0 [--build-only|--run-smoke]"
  exit 2
fi

echo "=== NPB-OMP smoke: CG Class S, OMP_NUM_THREADS=4 ==="
rm -rf "$ROOT/NPB3.4-OMP"/data_????????T??????
( cd "$ROOT/NPB3.4-OMP" && OMP_NUM_THREADS=4 OMP_PROC_BIND=true ./bin/cg.S.x ) | tee /tmp/npb_omp_smoke.out
grep -q 'Verification.*=.*SUCCESSFUL' /tmp/npb_omp_smoke.out
check_csv_dir "$ROOT/NPB3.4-OMP"

echo "=== NPB-MPI smoke: IS Class S, np=4 ==="
# Default suite conf uses PAPI_L2/L3_TCM (x86-oriented). On ARM those events are
# often unavailable; force a portable backend for the baseline IS smoke.
if [[ "$ARCH" == "aarch64" ]]; then
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" native perf_event_open "cpu-cycles,instructions" 2 0 0)
  trap 'restore_conf "'"$bak"'"; force_rebuild_measure "'"$ROOT/NPB3.4-MPI"'"; rm -f '"$SNAPSHOT" EXIT
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  (cd "$ROOT/NPB3.4-MPI" && make clean >/dev/null 2>&1 || true; make IS CLASS=S)
fi
rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/is.S.x ) | tee /tmp/npb_mpi_smoke.out
grep -q 'SUCCESSFUL' /tmp/npb_mpi_smoke.out
check_csv_dir "$ROOT/NPB3.4-MPI" 4
if [[ "$ARCH" == "aarch64" ]]; then
  restore_conf "$bak"
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  trap 'rm -f '"$SNAPSHOT" EXIT
fi

echo "=== NPB-MPI per-step: MG+CG Class S (ENABLE_PER_STEP_TIMING=1) ==="
if [[ "$ARCH" == "aarch64" ]]; then
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" native perf_event_open "cpu-cycles,instructions" 2 0 1)
else
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" native papi_read "PAPI_TOT_CYC,PAPI_TOT_INS" 2 0 1)
fi
trap 'restore_conf "'"$bak"'"; force_rebuild_measure "'"$ROOT/NPB3.4-MPI"'"; rm -f '"$SNAPSHOT" EXIT
force_rebuild_measure "$ROOT/NPB3.4-MPI"
(cd "$ROOT/NPB3.4-MPI" && make clean >/dev/null 2>&1 || true; make MG CLASS=S && make CG CLASS=S)
rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/mg.S.x ) | tee /tmp/npb_mpi_mg_step.out
grep -q 'VERIFICATION SUCCESSFUL\|Verification.*=.*SUCCESSFUL\|SUCCESSFUL' /tmp/npb_mpi_mg_step.out
check_csv_dir "$ROOT/NPB3.4-MPI" 4
check_per_step_mg_cg "$ROOT/NPB3.4-MPI"
rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/cg.S.x ) | tee /tmp/npb_mpi_cg_step.out
grep -q 'VERIFICATION SUCCESSFUL\|Verification.*=.*SUCCESSFUL\|SUCCESSFUL' /tmp/npb_mpi_cg_step.out
check_csv_dir "$ROOT/NPB3.4-MPI" 4
check_per_step_mg_cg "$ROOT/NPB3.4-MPI"

echo "=== NPB-MPI per-step OFF regression (ENABLE_PER_STEP_TIMING=0) ==="
restore_conf "$bak"
if [[ "$ARCH" == "aarch64" ]]; then
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" native perf_event_open "cpu-cycles,instructions" 2 0 0)
else
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" native papi_read "PAPI_TOT_CYC,PAPI_TOT_INS" 2 0 0)
fi
force_rebuild_measure "$ROOT/NPB3.4-MPI"
(cd "$ROOT/NPB3.4-MPI" && make clean >/dev/null 2>&1 || true; make MG CLASS=S)
rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/mg.S.x ) | tee /tmp/npb_mpi_mg_nostep.out
grep -q 'SUCCESSFUL' /tmp/npb_mpi_mg_nostep.out
dir=$(ls -dt "$ROOT/NPB3.4-MPI"/data_????????T?????? | head -1)
steps=$(awk -F, 'NR>1 && $6+0==1000 {n++} END{print n+0}' "$dir"/npb-mpi_mg_*.csv)
[[ "$steps" == "0" ]] || { echo "ERROR: expected 0 step rows with switch=0, got $steps"; exit 1; }
echo "per-step OFF OK"
restore_conf "$bak"
force_rebuild_measure "$ROOT/NPB3.4-MPI"
trap 'rm -f '"$SNAPSHOT" EXIT

# Platform-specific extra combo (rebuild with temp conf, then restore)
if [[ "$ARCH" == "x86_64" ]]; then
  echo "=== NPB-OMP combo: clock_gettime + perf_event_open ==="
  bak=$(apply_temp_conf "$ROOT/NPB3.4-OMP" clock_gettime perf_event_open "cpu-cycles,instructions" 2 0)
  trap 'restore_conf "'"$bak"'"; force_rebuild_measure "'"$ROOT/NPB3.4-OMP"'"; rm -f '"$SNAPSHOT" EXIT
  force_rebuild_measure "$ROOT/NPB3.4-OMP"
  (cd "$ROOT/NPB3.4-OMP" && make clean >/dev/null 2>&1 || true; make CG CLASS=S)
  rm -rf "$ROOT/NPB3.4-OMP"/data_????????T??????
  ( cd "$ROOT/NPB3.4-OMP" && OMP_NUM_THREADS=4 OMP_PROC_BIND=true ./bin/cg.S.x ) | tee /tmp/npb_omp_combo.out
  grep -q 'Verification.*=.*SUCCESSFUL' /tmp/npb_omp_combo.out
  check_csv_dir "$ROOT/NPB3.4-OMP"
  restore_conf "$bak"
  force_rebuild_measure "$ROOT/NPB3.4-OMP"
  trap 'rm -f '"$SNAPSHOT" EXIT

  echo "=== NPB-MPI combo: mpi_wtime + papi_read ==="
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" mpi_wtime papi_read "PAPI_TOT_CYC,PAPI_TOT_INS" 2 0)
  trap 'restore_conf "'"$bak"'"; force_rebuild_measure "'"$ROOT/NPB3.4-MPI"'"; rm -f '"$SNAPSHOT" EXIT
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  (cd "$ROOT/NPB3.4-MPI" && make clean >/dev/null 2>&1 || true; make IS CLASS=S)
  rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
  ( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/is.S.x ) | tee /tmp/npb_mpi_combo.out
  grep -q 'SUCCESSFUL' /tmp/npb_mpi_combo.out
  check_csv_dir "$ROOT/NPB3.4-MPI" 4
  restore_conf "$bak"
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  trap 'rm -f '"$SNAPSHOT" EXIT

  # restore native conf builds
  (cd "$ROOT/NPB3.4-OMP" && make CG CLASS=S)
  (cd "$ROOT/NPB3.4-MPI" && make IS CLASS=S)
elif [[ "$ARCH" == "aarch64" ]]; then
  NSTP="${TACVAR_NSTP_ARM:-10}"
  echo "=== NPB-OMP combo: cntvct_el0_dmb + asm (nstp=$NSTP) ==="
  bak=$(apply_temp_conf "$ROOT/NPB3.4-OMP" cntvct_el0_dmb asm "CPU_CYCLES,INST_RETIRED" 2 "$NSTP")
  trap 'restore_conf "'"$bak"'"; force_rebuild_measure "'"$ROOT/NPB3.4-OMP"'"; rm -f '"$SNAPSHOT" EXIT
  force_rebuild_measure "$ROOT/NPB3.4-OMP"
  (cd "$ROOT/NPB3.4-OMP" && make clean >/dev/null 2>&1 || true; make CG CLASS=S)
  rm -rf "$ROOT/NPB3.4-OMP"/data_????????T??????
  ( cd "$ROOT/NPB3.4-OMP" && OMP_NUM_THREADS=4 OMP_PROC_BIND=true ./bin/cg.S.x ) | tee /tmp/npb_omp_combo.out
  grep -q 'Verification.*=.*SUCCESSFUL' /tmp/npb_omp_combo.out
  check_csv_dir "$ROOT/NPB3.4-OMP"
  restore_conf "$bak"
  force_rebuild_measure "$ROOT/NPB3.4-OMP"

  echo "=== NPB-MPI combo: mpi_wtime + perf_event_open ==="
  # Note: on this HiSilicon host PAPI perf_event component is disabled
  # (libpfm4 no default PMU); use perf_event_open instead of papi_read.
  bak=$(apply_temp_conf "$ROOT/NPB3.4-MPI" mpi_wtime perf_event_open "cpu-cycles,instructions" 2 0)
  trap 'restore_conf "'"$bak"'"; force_rebuild_measure "'"$ROOT/NPB3.4-MPI"'"; rm -f '"$SNAPSHOT" EXIT
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  (cd "$ROOT/NPB3.4-MPI" && make clean >/dev/null 2>&1 || true; make IS CLASS=S)
  rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
  ( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/is.S.x ) | tee /tmp/npb_mpi_combo.out
  grep -q 'SUCCESSFUL' /tmp/npb_mpi_combo.out
  check_csv_dir "$ROOT/NPB3.4-MPI" 4
  restore_conf "$bak"
  force_rebuild_measure "$ROOT/NPB3.4-MPI"
  (cd "$ROOT/NPB3.4-OMP" && make CG CLASS=S)
  (cd "$ROOT/NPB3.4-MPI" && make IS CLASS=S)
fi

check_protected
echo "NPB smoke OK"
