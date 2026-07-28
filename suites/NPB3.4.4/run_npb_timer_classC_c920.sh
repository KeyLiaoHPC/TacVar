#!/usr/bin/env bash
# Compatibility wrapper for the legacy c920bn1 Class C campaign defaults.
#
# Prefer the two-stage workflow in utils/scripts/ (see utils/scripts/README.md and
# SKILLS/test-npb-mpi-timer-comparison/SKILL.md). This wrapper labels the
# historical c920/Class C defaults and either:
#   - delegates analysis to utils/scripts/analyze_npb_timer_campaign.py, or
#   - runs the preserved one-shot matrix via run_npb_timer_classC_c920.legacy.sh
#
# Legacy defaults (NOT silent defaults for new campaigns):
#   host=c920bn1, class=C, np=128 (BT/SP=121),
#   timers=native,clock_gettime,papi_get_real_nsec,cntvct_el0,cntvct_el0_dmb,
#   profiles=e0 (none) + e4 (papi_read, four PAPI events)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
TACVAR_ROOT="$(cd "$ROOT/../.." && pwd)"
SCRIPTS="$TACVAR_ROOT/utils/scripts"
LEGACY="$ROOT/run_npb_timer_classC_c920.legacy.sh"
PHASE="${1:-}"

case "$PHASE" in
  summarize|analyze)
    shift || true
    RESULTS="${1:-${RESUME_DIR:-}}"
    if [[ -z "$RESULTS" ]]; then
      RESULTS=$(ls -dt "$ROOT/NPB3.4-MPI"/results_c920bn1_classC_* 2>/dev/null | head -1 || true)
    else
      shift || true
    fi
    [[ -n "$RESULTS" && -d "$RESULTS" ]] || {
      echo "Usage: $0 summarize /path/to/results_..." >&2
      exit 2
    }
    PY=python3
    [[ -x "$SCRIPTS/.venv/bin/python" ]] && PY="$SCRIPTS/.venv/bin/python"
    echo "[wrapper] analyze via $SCRIPTS/analyze_npb_timer_campaign.py" >&2
    exec "$PY" "$SCRIPTS/analyze_npb_timer_campaign.py" "$RESULTS" --allow-partial "$@"
    ;;
  all|preflight|matrix|"")
    if [[ -x "$LEGACY" ]]; then
      echo "[wrapper] WARNING: invoking legacy one-shot Class C script." >&2
      echo "[wrapper] New campaigns should use utils/scripts/generate_npb_timer_jobs.py" >&2
      exec bash "$LEGACY" "$@"
    fi
    echo "Legacy script missing: $LEGACY" >&2
    exit 2
    ;;
  *)
    echo "Usage: $0 [all|preflight|matrix|summarize] ..." >&2
    exit 2
    ;;
esac
