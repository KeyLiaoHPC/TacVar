#!/usr/bin/env bash
# =============================================================================
# watch_dualtf.sh — read-only NFS progress for DualTF ROOT trees
#
#   bash scripts/watch_dualtf.sh [ROOT ...]
#   bash scripts/watch_dualtf.sh --loop 300 [ROOT ...]
#
# Default ROOT list: the three 20260815 trees under this workspace.
# Prints one line per host plus a SUM line. Does not SSH.
# =============================================================================
set -uo pipefail

INTERVAL=0
ROOTS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --loop)
      INTERVAL="${2:-300}"
      shift 2
      ;;
    --loop=*)
      INTERVAL="${1#--loop=}"
      shift
      ;;
    *)
      ROOTS+=("$1")
      shift
      ;;
  esac
done

WS="/mnt/keylabmain/nfs/hpckey/03-Project/TacVar_NPC"
if [[ ${#ROOTS[@]} -eq 0 ]]; then
  for h in c920bn1 cgnr6760pn1 camd9755n2; do
    ROOTS+=("$WS/TacVar_NPC_${h}/suites/NPB3.4.4/NPB3.4-MPI/matrix/${h}_dualtf_20260815")
  done
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

count_one() {
  local root="$1"
  local host
  host="$(basename "$root")"
  host="${host%%_dualtf_*}"

  local met_ok=0 tf_ok=0 fail=0 last="-"
  local combo kdir site_dir

  if [[ -d "$root" ]]; then
    for combo in "$root"/*; do
      [[ -d "$combo/met" ]] || continue
      [[ "$(basename "$combo")" == logs ]] && continue
      for kdir in "$combo"/met/*; do
        [[ -d "$kdir" ]] || continue
        [[ -f "$kdir/.met_ok" ]] && met_ok=$((met_ok + 1))
      done
      if [[ -d "$combo/tf" ]]; then
        for site_dir in "$combo"/tf/*; do
          [[ -d "$site_dir" ]] || continue
          [[ "$(basename "$site_dir")" == *_r*_l* ]] || continue
          if has_rank_csv "$site_dir/tfs" && has_rank_csv "$site_dir/tfe"; then
            tf_ok=$((tf_ok + 1))
          fi
        done
      fi
    done
    if [[ -d "$root/logs" ]]; then
      fail="$(grep -h -c ' FAIL ' "$root/logs"/*.log 2>/dev/null | awk '{s+=$1} END {print s+0}')"
      last="$(grep -h -E ' (OK|FAIL) ' "$root/logs"/*.log 2>/dev/null | tail -n 1 | tr -s ' ' | cut -c1-80)"
      [[ -z "$last" ]] && last="-"
    fi
  else
    last="ROOT missing"
  fi

  printf '%s met_ok=%d/36 (6x6 kernels) tf_ok=%d/240 (40x6 sites) fail=%d last=%s\n' \
    "$host" "$met_ok" "$tf_ok" "$fail" "$last"
  echo "$met_ok $tf_ok $fail"
}

print_report() {
  local tot_m=0 tot_t=0 tot_f=0
  local line rest m t f
  echo "[$(date '+%F %T')] DualTF NFS watch"
  for root in "${ROOTS[@]}"; do
    rest="$(count_one "$root")"
    line="$(printf '%s\n' "$rest" | head -n 1)"
    m="$(printf '%s\n' "$rest" | tail -n 1 | awk '{print $1}')"
    t="$(printf '%s\n' "$rest" | tail -n 1 | awk '{print $2}')"
    f="$(printf '%s\n' "$rest" | tail -n 1 | awk '{print $3}')"
    echo "  $line"
    tot_m=$((tot_m + m))
    tot_t=$((tot_t + t))
    tot_f=$((tot_f + f))
  done
  printf '  SUM met_ok=%d/108 tf_ok=%d/720 fail=%d\n' "$tot_m" "$tot_t" "$tot_f"
}

if [[ "$INTERVAL" -gt 0 ]]; then
  while true; do
    print_report
    sleep "$INTERVAL"
  done
else
  print_report
fi
