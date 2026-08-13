#!/usr/bin/env bash
# Full in-situ TF pipeline on KunPeng host c920bn1.
# Toolchain: ~/01-App/openmpi-5.0.8 and papi-7.2.0b2
# Launch: mpirun --map-by core --bind-to core -np 128
# Workload: CG CLASS=C, timer cntvct_el0, TACVAR_NSTP=10, counters none.
set -euo pipefail

HOST_EXPECT="${TACVAR_ARM_HOST:-c920bn1}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SUITE="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO="$(cd "$SUITE/../../.." && pwd)"
SELF="$SCRIPT_DIR/$(basename "$0")"

if [[ "$(hostname -s 2>/dev/null || hostname)" != "$HOST_EXPECT" ]]; then
  ssh -o BatchMode=yes -o ConnectTimeout=10 "$HOST_EXPECT" "test -d '$REPO'" \
    || { echo "ERROR: $HOST_EXPECT cannot see $REPO"; exit 1; }
  exec ssh "$HOST_EXPECT" "bash '$SELF'"
fi

# shellcheck disable=SC1091
source "$SCRIPT_DIR/setup_env.sh"

ARCH=$(uname -m)
[[ "$ARCH" == "aarch64" ]] || { echo "ERROR: expected aarch64, got $ARCH"; exit 1; }

NP="${NP:-128}"
BENCH=cg
CLASS=C
CONF="$SUITE/tacvar.conf"
BAK="$SUITE/tacvar.conf.c920bn1.bak"
PROTECT="$REPO/src/measure/tests/protected_sources.txt"
SNAPSHOT=$(mktemp)

snapshot_protected() {
  (cd "$REPO" && while IFS= read -r p; do
    [[ -z "$p" || "$p" =~ ^# ]] && continue
    if [[ -f "$p" ]]; then sha256sum "$p"
    elif [[ -d "$p" ]]; then
      find "$p" -type f \( -name '*.f' -o -name '*.f90' -o -name '*.c' -o -name '*.h' \) \
        ! -name 'npbparams.h' ! -name 'mpinpb.f90' ! -name 'mpinpb.h' -print0 \
        | sort -z | xargs -0 sha256sum 2>/dev/null || true
    fi
  done < "$PROTECT")
}

restore_conf() {
  if [[ -f "$BAK" ]]; then
    mv -f "$BAK" "$CONF"
    touch "$CONF"
  fi
}

trap 'restore_conf; rm -f "$SNAPSHOT"' EXIT
snapshot_protected >"$SNAPSHOT"
cp "$CONF" "$BAK"

cat >"$CONF" <<EOF
TACVAR_TIMER=cntvct_el0
TACVAR_NSTP=${TACVAR_NSTP_ARM:-10}
TACVAR_COUNTER_BACKEND=none
TACVAR_COUNTER_COUNT=0
TACVAR_COUNTER_NAMES=
TACVAR_OUTPUT_ROOT=.
TACVAR_TF_SAMPLING_MODE=OFF
TACVAR_TF_DATA_ROOT=
TACVAR_TF_NSPG=0
EOF
touch "$CONF"

export NP BENCH CLASS
"$SCRIPT_DIR/run_npb_measure.sh" all

DATA=$(ls -dt "$SUITE"/data_????????T?????? 2>/dev/null | head -1)
[[ -n "$DATA" ]] || { echo "ERROR: no data_*"; exit 1; }
[[ -f "$DATA/median.csv" ]] || { echo "ERROR: no median.csv"; exit 1; }
grep -q '^CG,C,' "$DATA/median.csv" || { echo "ERROR: median.csv missing CG,C rows"; exit 1; }

ncsv() { find "$1" -maxdepth 1 -name '*_r[0-9]*_t[0-9]*_p[0-9]*.csv' | wc -l; }
[[ -d "$DATA/cg.C" ]] || { echo "ERROR: missing $DATA/cg.C"; exit 1; }
[[ "$(ncsv "$DATA/cg.C")" -eq "$NP" ]] || {
  echo "ERROR: cg.C rank CSV count $(ncsv "$DATA/cg.C") != $NP"; exit 1
}
[[ -x "$SUITE/bin/cg.C_tf.x" ]] || {
  echo "ERROR: missing TF binary cg.C_tf.x"; exit 1
}
[[ "$(ncsv "$DATA/cg.C_tf")" -eq "$NP" ]] || {
  echo "ERROR: cg.C_tf rank CSV count $(ncsv "$DATA/cg.C_tf") != $NP"; exit 1
}
shopt -s nullglob
filt_dirs=("$DATA"/cg.C_filt/r*_l*)
[[ ${#filt_dirs[@]} -ge 1 ]] || { echo "ERROR: no cg.C_filt/r*_l* dirs"; exit 1; }
head -1 "$DATA/median.csv" | grep -q ',ngauge' || {
  echo "ERROR: median.csv missing ngauge column"; exit 1
}
tfcdf=$(find "$DATA/cg.C_filt" -name tf_cdf.csv 2>/dev/null | head -1)
[[ -n "$tfcdf" ]] || { echo "ERROR: no tf_cdf.csv under cg.C_filt"; exit 1; }
nlines=$(wc -l < "$tfcdf")
[[ "$nlines" -eq 1001 ]] || {
  echo "ERROR: $tfcdf has $nlines lines, expected 1001"; exit 1
}

now=$(mktemp)
snapshot_protected >"$now"
if ! diff -q "$SNAPSHOT" "$now" >/dev/null; then
  echo "ERROR: protected NPB sources were modified"
  diff -u "$SNAPSHOT" "$now" || true
  rm -f "$now"
  exit 1
fi
rm -f "$now"
echo "c920bn1 TF test OK  data=$DATA"
