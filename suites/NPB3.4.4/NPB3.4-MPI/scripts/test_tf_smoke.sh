#!/usr/bin/env bash
# Local/CI TF smoke: CG CLASS=S, np=2, native timer, no counters.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SUITE="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO="$(cd "$SUITE/../../.." && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/setup_env.sh"

NP="${NP:-2}"
BENCH=cg
CLASS=S
CONF="$SUITE/tacvar.conf"
BAK="$SUITE/tacvar.conf.tfsmoke.bak"
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

NSTP=0
TIMER=native
if [[ "$(uname -m)" == "aarch64" ]]; then
  TIMER=cntvct_el0
  NSTP="${TACVAR_NSTP_ARM:-10}"
fi

cat >"$CONF" <<EOF
TACVAR_TIMER=$TIMER
TACVAR_NSTP=$NSTP
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
[[ -x "$SUITE/bin/cg.S_tf.x" ]] || {
  echo "ERROR: missing cg.S_tf.x"; exit 1
}
[[ -d "$DATA/cg.S" && -d "$DATA/cg.S_tf" ]] || {
  echo "ERROR: missing cg.S / cg.S_tf dirs under $DATA"; exit 1
}
head -1 "$DATA/median.csv" | grep -q ',ngauge' || {
  echo "ERROR: median.csv missing ngauge column"; exit 1
}
tfcdf=$(find "$DATA/cg.S_filt" -name tf_cdf.csv 2>/dev/null | head -1)
[[ -n "$tfcdf" ]] || { echo "ERROR: no tf_cdf.csv under cg.S_filt"; exit 1; }
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
echo "TF smoke OK  data=$DATA"
