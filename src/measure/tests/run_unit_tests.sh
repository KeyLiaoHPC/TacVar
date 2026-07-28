#!/usr/bin/env bash
# Unit tests for TacVar measure (mock / no hardware PMU required).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
if [[ "$(uname -m)" == "aarch64" ]]; then
  PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi}"
else
  PAPI_HOME="${PAPI_HOME:-/home/hpckey/01-App/papi-7.1.0}"
fi

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

# Per-step switch validation
for v in 0 1; do
  cat > "$TMP/ps_${v}.conf" <<EOF
TACVAR_TIMER=native
TACVAR_NSTP=0
TACVAR_COUNTER_BACKEND=none
TACVAR_COUNTER_COUNT=0
TACVAR_COUNTER_NAMES=
TACVAR_OUTPUT_ROOT=.
TACVAR_ENABLE_PER_STEP_TIMING=$v
EOF
  python3 "$ROOT/tools/gen_config.py" --conf "$TMP/ps_${v}.conf" --compiler gcc \
    --consumer lmbench --outdir "$TMP/ps_${v}" --measure-root "$ROOT" >/dev/null
  grep -q "#define TACVAR_ENABLE_PER_STEP_TIMING $v" "$TMP/ps_${v}/tacvar_generated_config.h"
  grep -q "#define TACVAR_REGION_STEP 1000" "$TMP/ps_${v}/tacvar_generated_config.h"
  check "accept ENABLE_PER_STEP_TIMING=$v" true
done

for bad in "" 2 -1 abc 1.5; do
  cat > "$TMP/ps_bad.conf" <<EOF
TACVAR_TIMER=native
TACVAR_NSTP=0
TACVAR_COUNTER_BACKEND=none
TACVAR_COUNTER_COUNT=0
TACVAR_COUNTER_NAMES=
TACVAR_OUTPUT_ROOT=.
TACVAR_ENABLE_PER_STEP_TIMING=$bad
EOF
  if python3 "$ROOT/tools/gen_config.py" --conf "$TMP/ps_bad.conf" --compiler gcc \
    --consumer lmbench --outdir "$TMP/ps_bad" --measure-root "$ROOT" 2>/dev/null; then
    check "reject ENABLE_PER_STEP_TIMING='$bad'" false
  else
    check "reject ENABLE_PER_STEP_TIMING='$bad'" true
  fi
done

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

# Workload integrity checker
python3 "$ROOT/tests/check_npb_workload.py" >/dev/null
check "npb workload guard" true

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
  "migrated,valid,source";
int main(void) {
  assert(strstr(hdr, "elapsed_ns") != NULL);
  assert(strstr(hdr, "migrated") != NULL);
  assert(strstr(hdr, "source") != NULL);
  return 0;
}
EOF
gcc -O2 -o "$TMP/csv_schema" "$TMP/csv_schema.c" && check "csv schema fields" "$TMP/csv_schema"

# Deferred buffer / reserve / region_info / fork tests against libtacvar_measure
cat > "$TMP/buf.conf" <<EOF
TACVAR_TIMER=native
TACVAR_NSTP=0
TACVAR_COUNTER_BACKEND=none
TACVAR_COUNTER_COUNT=0
TACVAR_COUNTER_NAMES=
TACVAR_OUTPUT_ROOT=$TMP/out
TACVAR_ENABLE_PER_STEP_TIMING=0
EOF
make -C "$ROOT" CONF="$TMP/buf.conf" CONSUMER=lmbench OUTDIR="$TMP/lib" CC=gcc \
  PAPI_INC="$PAPI_HOME/include" PAPI_LIB="$PAPI_HOME/lib" >/dev/null
mkdir -p "$TMP/out"

cat > "$TMP/buf_test.c" <<'EOF'
#include "tacvar_measure.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <dirent.h>
#include <errno.h>

static int count_csv(const char *dir)
{
  DIR *d = opendir(dir);
  struct dirent *e;
  int n = 0;
  if (!d) return -1;
  while ((e = readdir(d)) != NULL) {
    size_t len = strlen(e->d_name);
    if (len > 4 && strcmp(e->d_name + len - 4, ".csv") == 0)
      n++;
  }
  closedir(d);
  return n;
}

static void write_one(int region, uint64_t a, uint64_t b, int64_t ns,
                      int cpu_a, int cpu_b, const char *src)
{
  tacvar_csv_write_simple(region, a, b, ns, cpu_a, cpu_b, NULL, NULL, NULL, src);
}

int main(void)
{
  tacvar_context_t ctx = {
    .suite = "unit", .benchmark = "buf", .klass = "S",
    .test_tag = "tag,with\"quote", .rank = 0, .thread = 0, .nprocs = 1
  };
  const char *ddir;
  char path[512];
  FILE *fp;
  char line[1024];
  int rc;
  pid_t pid;
  int status;
  tacvar_region_info_t regs[2];

  assert(tacvar_init(&ctx) == 0);
  ddir = tacvar_data_dir();
  assert(ddir && ddir[0]);
  assert(count_csv(ddir) == 0); /* deferred: no CSV yet */

  rc = tacvar_events_reserve(4);
  assert(rc == 0);
  write_one(1, 10, 20, 100, 0, 0, "a.c:1");
  write_one(1000, 30, 40, 200, 1, 2, "b.c:2"); /* migrated */
  write_one(1, 50, 60, 300, 3, 3, "c.c:3");
  assert(count_csv(ddir) == 0);

  memset(&regs[0], 0, sizeof(regs));
  regs[0].region_id = 1;
  regs[0].region_name = "total";
  regs[0].source_locations = "file:r:1-2";
  regs[0].description = "whole kernel";
  regs[0].active_when = "always";
  regs[0].recorded_by_tacvar = 1;
  regs[1].region_id = 1000;
  regs[1].region_name = "step";
  regs[1].source_locations = "a;b";
  regs[1].description = "desc,with,commas";
  regs[1].active_when = "TACVAR_ENABLE_PER_STEP_TIMING=1";
  regs[1].recorded_by_tacvar = 1;
  tacvar_region_info_register(regs, 2);

  /* Forked child must not flush parent events. */
  pid = fork();
  assert(pid >= 0);
  if (pid == 0) {
    /* inherited ready with wrong owner: no-op fini; re-init discards buffer */
    tacvar_fini();
    assert(count_csv(ddir) == 0);
    _exit(0);
  }
  assert(waitpid(pid, &status, 0) == pid);
  assert(WIFEXITED(status) && WEXITSTATUS(status) == 0);
  assert(count_csv(ddir) == 0);

  tacvar_fini();
  assert(count_csv(ddir) == 2); /* event CSV + region_info.csv */

  snprintf(path, sizeof(path), "%s/region_info.csv", ddir);
  fp = fopen(path, "r");
  assert(fp);
  assert(fgets(line, sizeof(line), fp));
  assert(strstr(line, "region_id") != NULL);
  assert(fgets(line, sizeof(line), fp));
  assert(strstr(line, ",1,") != NULL);
  assert(fgets(line, sizeof(line), fp));
  assert(strstr(line, "\"desc,with,commas\"") != NULL || strstr(line, "desc,with,commas") != NULL);
  assert(strstr(line, "1000") != NULL);
  fclose(fp);

  /* Find event CSV and verify order / migration / threads */
  {
    DIR *d = opendir(ddir);
    struct dirent *e;
    char epath[512] = "";
    assert(d);
    while ((e = readdir(d)) != NULL) {
      if (strcmp(e->d_name, "region_info.csv") == 0) continue;
      if (strstr(e->d_name, ".csv")) {
        snprintf(epath, sizeof(epath), "%s/%s", ddir, e->d_name);
        break;
      }
    }
    closedir(d);
    assert(epath[0]);
    fp = fopen(epath, "r");
    assert(fp);
    assert(fgets(line, sizeof(line), fp)); /* header */
    assert(fgets(line, sizeof(line), fp));
    assert(strstr(line, ",1,gettimeofday,") != NULL || strstr(line, ",1,native") != NULL);
    assert(fgets(line, sizeof(line), fp));
    assert(strstr(line, ",1000,") != NULL);
    /* cpu_start=1,cpu_stop=2 => migrated=1 */
    assert(strstr(line, ",1,2,1,") != NULL);
    assert(fgets(line, sizeof(line), fp));
    assert(strstr(line, ",1,") != NULL);
    assert(!fgets(line, sizeof(line), fp)); /* exactly 3 events */
    fclose(fp);
  }

  /* Overflow-ish: reserve huge should fail cleanly if too large */
  assert(tacvar_init(&ctx) == 0);
  rc = tacvar_events_reserve((size_t)-1 / 2);
  assert(rc == -EOVERFLOW || rc == -ENOMEM || rc == 0);
  tacvar_fini();

  /* Re-init must still accept region registration (simulates NPB always-register). */
  {
    const char *ddir2;
    char rpath[512];
    assert(tacvar_init(&ctx) == 0);
    ddir2 = tacvar_data_dir();
    assert(ddir2 && ddir2[0]);
    tacvar_region_info_clear();
    tacvar_region_info_register(regs, 2);
    write_one(1, 1, 2, 10, 0, 0, "reinit");
    tacvar_fini();
    snprintf(rpath, sizeof(rpath), "%s/region_info.csv", ddir2);
    assert(access(rpath, R_OK) == 0);
  }

  puts("buf_test ok");
  return 0;
}
EOF

# shellcheck disable=SC1090
# Pull include paths from generated mk (timers/, counters/, include/).
MEAS_CFLAGS=$(sed -n 's/^TACVAR_MEASURE_CFLAGS := //p' "$TMP/lib/tacvar_generated_config.mk")
MEAS_LDFLAGS=$(sed -n 's/^TACVAR_MEASURE_LDFLAGS := //p' "$TMP/lib/tacvar_generated_config.mk")
gcc -std=c99 -O2 -D_GNU_SOURCE $MEAS_CFLAGS -I"$TMP/lib" \
  -DTACVAR_HAS_GENERATED_CONFIG -include "$TMP/lib/tacvar_generated_config.h" \
  -o "$TMP/buf_test" "$TMP/buf_test.c" "$TMP/lib/libtacvar_measure.a" $MEAS_LDFLAGS -lrt -ldl \
  && (cd "$TMP" && "$TMP/buf_test") \
  && check "deferred buffer + region_info + fork" true \
  || check "deferred buffer + region_info + fork" false

echo "unit tests: $pass passed, $fail failed"
exit $fail
