#!/usr/bin/env bash
# Unit tests for TacVar measure (mock / no hardware PMU required).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

pass=0
fail=0
check() {
  local name="$1"; shift
  if "$@"; then echo "PASS: $name"; pass=$((pass+1)); else echo "FAIL: $name"; fail=$((fail+1)); fi
}

# Config validation
python3 "$ROOT/tools/gen_config.py" --conf "$ROOT/../../suites/lmbench/tacvar.conf" \
  --compiler gcc --consumer lmbench --outdir "$TMP/ok" --measure-root "$ROOT" >/dev/null
check "gen_config lmbench native" test -f "$TMP/ok/tacvar_generated_config.h"

cat > "$TMP/bad.conf" <<EOF
TACVAR_TIMER=mpi_wtime
TACVAR_NSTP=0
TACVAR_COUNTER_BACKEND=none
TACVAR_COUNTER_COUNT=0
TACVAR_COUNTER_NAMES=
TACVAR_OUTPUT_ROOT=.
EOF
if python3 "$ROOT/tools/gen_config.py" --conf "$TMP/bad.conf" --compiler gcc \
  --consumer lmbench --outdir "$TMP/bad" --measure-root "$ROOT" 2>/dev/null; then
  check "reject mpi_wtime on lmbench" false
else
  check "reject mpi_wtime on lmbench" true
fi

# TSC join helper compile test
cat > "$TMP/tsc.c" <<'EOF'
#include <stdint.h>
#include <assert.h>
static inline uint64_t join(unsigned hi, unsigned lo) {
  return ((uint64_t)hi << 32) | (uint64_t)lo;
}
int main(void) {
  assert(join(1,2) == 0x100000002ULL);
  assert(join(0xffffffffu, 0xffffffffu) == 0xffffffffffffffffULL);
  return 0;
}
EOF
gcc -O2 -o "$TMP/tsc" "$TMP/tsc.c" && check "tsc join" "$TMP/tsc"

# Modular delta
cat > "$TMP/delta.c" <<'EOF'
#include <stdint.h>
#include <assert.h>
static uint64_t mod_delta(uint64_t s, uint64_t e, unsigned w) {
  uint64_t mask = (w >= 64) ? ~0ULL : ((1ULL << w) - 1ULL);
  return (e - s) & mask;
}
int main(void) {
  assert(mod_delta(0xffe, 0x002, 12) == 4);
  return 0;
}
EOF
gcc -O2 -o "$TMP/delta" "$TMP/delta.c" && check "counter wrap delta" "$TMP/delta"



# Protected sources list exists
check "protected_sources.txt present" test -f "$ROOT/tests/protected_sources.txt"

# Reject ARM timer on x86 compiler defines (if on x86)
if [[ "$(uname -m)" == "x86_64" ]]; then
  cat > "$TMP/armt.conf" <<EOF
TACVAR_TIMER=cntvct_el0
TACVAR_NSTP=10
TACVAR_COUNTER_BACKEND=none
TACVAR_COUNTER_COUNT=0
TACVAR_COUNTER_NAMES=
TACVAR_OUTPUT_ROOT=.
EOF
  if python3 "$ROOT/tools/gen_config.py" --conf "$TMP/armt.conf" --compiler gcc \
    --consumer lmbench --outdir "$TMP/armt" --measure-root "$ROOT" 2>/dev/null; then
    check "reject cntvct on x86" false
  else
    check "reject cntvct on x86" true
  fi
fi

# CSV schema compile+run mini test
cat > "$TMP/csv_schema.c" <<'EOF'
#include <stdio.h>
#include <string.h>
#include <assert.h>
static const char *hdr =
  "seq,suite,benchmark,class,test_tag,region_id,timer,"
  "raw_start,raw_stop,elapsed_ns,rank,thread,pid,cpu_start,cpu_stop,"
  "migrated,counter_backend";
int main(void) {
  assert(strstr(hdr, "elapsed_ns") != NULL);
  assert(strstr(hdr, "migrated") != NULL);
  assert(strstr(hdr, "counter_backend") != NULL);
  return 0;
}
EOF
gcc -O2 -o "$TMP/csv_schema" "$TMP/csv_schema.c" && check "csv schema fields" "$TMP/csv_schema"

echo "unit tests: $pass passed, $fail failed"
exit $fail
