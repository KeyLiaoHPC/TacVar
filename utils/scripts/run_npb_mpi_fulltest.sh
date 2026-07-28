#!/usr/bin/env bash
# Run a generated NPB-MPI campaign fulltest (01_fulltest.sh) after pretest PASS.
# Usage:
#   run_npb_mpi_fulltest.sh /path/to/results_...
#   run_npb_mpi_fulltest.sh --results-dir /path/to/results_... [--generate]
set -euo pipefail

SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
RESULTS=""
GENERATE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS="$2"; shift 2 ;;
    --generate) GENERATE=1; shift ;;
    -h|--help)
      sed -n '1,6p' "$0"
      exit 0
      ;;
    *)
      if [[ -z "$RESULTS" && -d "$1" ]]; then
        RESULTS="$1"; shift
      else
        echo "Unknown arg: $1" >&2; exit 2
      fi
      ;;
  esac
done

[[ -n "$RESULTS" ]] || { echo "Usage: $0 /path/to/results_..." >&2; exit 2; }

if [[ "$GENERATE" -eq 1 || ! -x "$RESULTS/01_fulltest.sh" ]]; then
  PY=python3
  [[ -x "$SCRIPTS/.venv/bin/python" ]] && PY="$SCRIPTS/.venv/bin/python"
  "$PY" "$SCRIPTS/generate_npb_timer_jobs.py" fulltest --results-dir "$RESULTS"
fi

[[ -x "$RESULTS/01_fulltest.sh" ]] || { echo "Missing $RESULTS/01_fulltest.sh" >&2; exit 2; }
exec bash "$RESULTS/01_fulltest.sh"
