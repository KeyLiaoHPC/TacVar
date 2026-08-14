#!/usr/bin/env bash
# Build-smoke for NPB-MPI in-situ TF binaries (does not edit protected kernels).
# Usage: scripts/test_tf_insitu.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SUITE="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/setup_env.sh"
CONF="$SUITE/tacvar.conf"
BAK=$(mktemp)
DATA=$(mktemp -d)
trap 'mv -f "$BAK" "$CONF"; rm -rf "$DATA"' EXIT

cp "$CONF" "$BAK"

# Dummy measurement products so TF init can load ngauge.
mkdir -p "$DATA"
echo "1.0" >"$DATA/nspg.txt"
cat >"$DATA/median.csv" <<EOF
kernel,class,region_id,loc_id,median,ngauge
IS,S,0,1,10000,10000
EOF

python3 - <<PY
from pathlib import Path
p = Path("$CONF")
text = p.read_text()
keys = {
    "TACVAR_TF_SAMPLING_MODE": "ON",
    "TACVAR_TF_DATA_ROOT": "$DATA",
    "TACVAR_TF_REG_ID": "0",
    "TACVAR_TF_LOC_ID": "1",
}
lines = []
seen = set()
for line in text.splitlines():
    if "=" in line and not line.strip().startswith("#"):
        k = line.split("=", 1)[0].strip()
        if k in keys:
            lines.append(f"{k}={keys[k]}")
            seen.add(k)
            continue
    lines.append(line)
for k, v in keys.items():
    if k not in seen:
        lines.append(f"{k}={v}")
p.write_text("\n".join(lines) + "\n")
PY

make -C "$SUITE" tacvar_clean
make -C "$SUITE" IS CLASS=S

test -x "$SUITE/bin/is.S_tfs.x"
test -x "$SUITE/bin/is.S_tfe.x"
echo "PASS: $SUITE/bin/is.S_tfs.x and is.S_tfe.x"

# Dummy FilT inputs: original met + tfs/tfe CSVs under a separate tree.
FILT_DATA=$(mktemp -d)
trap 'mv -f "$BAK" "$CONF"; rm -rf "$DATA" "$FILT_DATA"' EXIT
hdr="seq,region_id,loc_id,raw_start,raw_stop,elapsed_ns,rank,thread,pid,cpu_start,cpu_stop,migrated,valid"
write_csv() {
    local dir="$1" elapsed="$2"
    mkdir -p "$dir"
    {
        echo "$hdr"
        i=1
        while [ "$i" -le 40 ]; do
            echo "$i,0,1,0,$elapsed,$elapsed,0,0,1,0,0,0,1"
            i=$((i + 1))
        done
    } >"$dir/rank0.csv"
}
echo "1.0" >"$FILT_DATA/nspg.txt"
write_csv "$FILT_DATA/is.S" 5000
write_csv "$FILT_DATA/is.S_tfs" 6200
write_csv "$FILT_DATA/is.S_tfe" 6100
python3 "$SCRIPT_DIR/get_met_stat.py" "$FILT_DATA"
grep -q 'IS,S,0,1,5000,5000' "$FILT_DATA/median.csv" || {
    echo "FAIL: ngauge != round(5000/1.0)"
    cat "$FILT_DATA/median.csv"
    exit 1
}
make -C "$SUITE" filt
python3 "$SCRIPT_DIR/run_filt.py" "$FILT_DATA" --filt "$SUITE/bin/filt.x" --nsamp 2000 --width 50
test -f "$FILT_DATA/is.S_filt/r0_l1/met.csv"
test -f "$FILT_DATA/is.S_filt/r0_l1/tf.csv"
test -f "$FILT_DATA/is.S_filt/r0_l1/tr_hist.csv" || test -f "$FILT_DATA/is.S_filt/r0_l1/er.out"
echo "PASS: get_met_stat.py + run_filt.py -> $FILT_DATA/is.S_filt/r0_l1"
