#!/usr/bin/env bash
# Run a generated NPB-MPI campaign pretest (00_pretest.sh) from a results directory.
# Usage:
#   run_npb_mpi_pretest.sh /path/to/results_...
#   run_npb_mpi_pretest.sh --spec /path/to/campaign.json [--results-dir DIR]
set -euo pipefail

SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
SPEC=""
RESULTS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --spec) SPEC="$2"; shift 2 ;;
    --results-dir) RESULTS="$2"; shift 2 ;;
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

if [[ -n "$SPEC" ]]; then
  PY=python3
  [[ -x "$SCRIPTS/.venv/bin/python" ]] && PY="$SCRIPTS/.venv/bin/python"
  args=(pretest --spec "$SPEC")
  [[ -n "$RESULTS" ]] && args+=(--results-dir "$RESULTS")
  mapfile -t out < <("$PY" "$SCRIPTS/generate_npb_timer_jobs.py" "${args[@]}")
  RESULTS="${out[0]}"
fi

[[ -n "$RESULTS" ]] || { echo "Usage: $0 /path/to/results_... | --spec campaign.json" >&2; exit 2; }
[[ -x "$RESULTS/00_pretest.sh" ]] || { echo "Missing $RESULTS/00_pretest.sh — generate first" >&2; exit 2; }
exec bash "$RESULTS/00_pretest.sh"
