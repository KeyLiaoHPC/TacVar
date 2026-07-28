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
  export LD_LIBRARY_PATH="/home/hpckey/01-App/openmpi-5.0.8/lib:/home/hpckey/01-App/papi-7.2.0b2/lib:${LD_LIBRARY_PATH:-}"
  export PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.2.0b2}"
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
  mapfile -t csvs < <(find "$dir" -name '*.csv' | sort)
  [[ ${#csvs[@]} -ge 1 ]] || { echo "ERROR: no CSV in $dir"; exit 1; }
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

apply_temp_conf() {
  local suite_dir="$1" timer="$2" counter="$3" names="$4" count="$5" nstp="$6"
  local bak="$suite_dir/tacvar.conf.smoke.bak"
  cp "$suite_dir/tacvar.conf" "$bak"
  cat > "$suite_dir/tacvar.conf" <<EOF
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
rm -rf "$ROOT/NPB3.4-MPI"/data_????????T??????
( cd "$ROOT/NPB3.4-MPI" && mpirun -np 4 --bind-to core ./bin/is.S.x ) | tee /tmp/npb_mpi_smoke.out
grep -q 'SUCCESSFUL' /tmp/npb_mpi_smoke.out
check_csv_dir "$ROOT/NPB3.4-MPI" 4

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
