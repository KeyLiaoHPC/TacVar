#!/usr/bin/env bash
# =============================================================================
# check_dualtf_complete.sh — gate: nspt/nspg + 36 met_ok + 240 isolated tfs/tfe
#
#   bash scripts/check_dualtf_complete.sh [ROOT]
# ROOT empty: source dualtf_common.sh and resolve_root.
# Exit 0 only when the tree is complete enough to filt.
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ROOT_ARG="${1:-}"
if [[ -n "$ROOT_ARG" ]]; then
  ROOT="$ROOT_ARG"
else
  # shellcheck source=dualtf_common.sh
  source "$SCRIPT_DIR/dualtf_common.sh"
  resolve_root || exit 2
fi

has_rank_csv() {
  local d="$1" f
  [[ -d "$d" ]] || return 1
  for f in "$d"/*.csv; do
    [[ -e "$f" ]] || continue
    [[ "$(basename "$f")" == "timer_info.csv" ]] && continue
    return 0
  done
  return 1
}

if [[ ! -f "$ROOT/nspt.txt" || ! -f "$ROOT/nspg.txt" ]]; then
  echo "incomplete: missing nspt.txt or nspg.txt under $ROOT" >&2
  exit 1
fi

met_ok=0
tf_ok=0
combos=0
for combo in "$ROOT"/*; do
  [[ -d "$combo/met" && -d "$combo/tf" ]] || continue
  [[ "$(basename "$combo")" == logs ]] && continue
  combos=$((combos + 1))
  for kdir in "$combo"/met/*; do
    [[ -d "$kdir" ]] || continue
    [[ -f "$kdir/.met_ok" ]] && met_ok=$((met_ok + 1))
  done
  for site_dir in "$combo"/tf/*; do
    [[ -d "$site_dir" ]] || continue
    [[ "$(basename "$site_dir")" == *_r*_l* ]] || continue
    if has_rank_csv "$site_dir/tfs" && has_rank_csv "$site_dir/tfe"; then
      tf_ok=$((tf_ok + 1))
    fi
  done
done

echo "ROOT=$ROOT combos=$combos met_ok=$met_ok/36 tf_ok=$tf_ok/240"
if [[ "$met_ok" -ge 36 && "$tf_ok" -ge 240 ]]; then
  echo "complete"
  exit 0
fi
echo "incomplete" >&2
exit 1
