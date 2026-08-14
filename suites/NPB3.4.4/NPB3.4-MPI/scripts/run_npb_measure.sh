#!/usr/bin/env bash
# Build and run NPB-MPI normal timing, nspg fit, and median.csv.
# Usage: scripts/run_npb_measure.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SUITE="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/setup_env.sh"

BENCH="${BENCH:-cg}"
CLASS="${CLASS:-S}"
NP="${NP:-2}"
MPIRUN_BIN="${MPIRUN_BIN:-mpirun}"
MPIRUN_ARGS=(--map-by core --bind-to core -np "$NP")
CONF="$SUITE/tacvar.conf"

BENCH_U=$(echo "$BENCH" | tr '[:lower:]' '[:upper:]')
BIN_BASE="$SUITE/bin/${BENCH}.${CLASS}"

set_conf_kv() {
  local key="$1" val="$2"
  if grep -q "^${key}=" "$CONF"; then
    sed -i "s|^${key}=.*|${key}=${val}|" "$CONF"
  else
    echo "${key}=${val}" >>"$CONF"
  fi
}

latest_data_dir() {
  ls -dt "$SUITE"/data_????????T?????? 2>/dev/null | head -1
}

build_nspg() {
  make -C "$SUITE" nspg PAPI_HOME="${PAPI_HOME}" NSPG_MPICC="$(command -v mpicc)"
}

run_nspg() {
  local dest="${1:-}"
  local tmp
  tmp=$(mktemp)
  (cd "$SUITE" && "$MPIRUN_BIN" "${MPIRUN_ARGS[@]}" "$SUITE/bin/test_nspg.x" "$tmp")
  if [[ -n "$dest" ]]; then
    mkdir -p "$(dirname "$dest")"
    cp "$tmp" "$dest"
  else
    cp "$tmp" "$SUITE/nspg.txt"
  fi
  rm -f "$tmp"
}

set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
set_conf_kv TACVAR_TF_DATA_ROOT ""
(cd "$SUITE" && make tacvar_clean && make "$BENCH_U" CLASS="$CLASS")
build_nspg
rm -rf "$SUITE"/data_????????T??????
(cd "$SUITE" && NPB_TIMER_FLAG=1 "$MPIRUN_BIN" "${MPIRUN_ARGS[@]}" "${BIN_BASE}.x" | tee /tmp/npb_measure.out)
grep -E 'Verification[[:space:]]*=[[:space:]]*SUCCESSFUL|SUCCESSFUL' /tmp/npb_measure.out >/dev/null \
  || { echo "NPB verification failed"; exit 1; }
d=$(latest_data_dir)
[[ -n "$d" ]] || { echo "no data_* after measure"; exit 1; }
run_nspg "$d/nspg.txt"
python3 "$SCRIPT_DIR/get_median.py" "$d"
echo "DATA_DIR=$d"
