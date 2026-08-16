#!/usr/bin/env bash
# =============================================================================
# run_dualtf_filt_workers.sh — parallel FilT, one exclusive claim per site
#
#   bash scripts/run_dualtf_filt_workers.sh [NWORKERS]
#
# Two workers never run filt.x on the same (combo, kernel, rid, lid).
# Skips sites that already have filter/<site>/wd.out or an existing claim.
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dualtf_common.sh
source "$SCRIPT_DIR/dualtf_common.sh"

resolve_root || exit 2

NWORKERS="${1:-8}"
if [[ ! "$NWORKERS" =~ ^[1-9][0-9]*$ ]]; then
  log "NWORKERS must be a positive integer"
  exit 2
fi
if [[ ! -x "$SUITE/bin/filt.x" ]]; then
  log "missing $SUITE/bin/filt.x; run make filt on cicx6330n1"
  exit 2
fi

CLAIM="$ROOT/logs/filt_claim"
LOG="$ROOT/logs/filt_all.log"
mkdir -p "$CLAIM" "$ROOT/logs"

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

claim_name() {
  local combo="$1" site="$2"
  local s="${combo}__${site}"
  echo "${s//+/_}"
}

# Claim any site already being processed under this ROOT.
claim_inflight() {
  local line combo_dir combo rid lid kern site sdir
  while IFS= read -r line; do
    [[ "$line" == *"$ROOT/"* ]] || continue
    combo_dir="$(sed -n 's/.*run_filt.py \([^ ]*\) .*/\1/p' <<<"$line")"
    [[ -n "$combo_dir" ]] || continue
    combo="$(basename "$combo_dir")"
    rid="$(sed -n 's/.*--rid \([0-9][0-9]*\).*/\1/p' <<<"$line")"
    lid="$(sed -n 's/.*--lid \([0-9][0-9]*\).*/\1/p' <<<"$line")"
    kern="$(sed -n 's/.*--kernel \([^ ]*\).*/\1/p' <<<"$line")"
    [[ -n "$rid" && -n "$lid" ]] || continue
    if [[ -n "$kern" ]]; then
      site="${kern}.C_r${rid}_l${lid}"
      mkdir -p "$CLAIM/$(claim_name "$combo" "$site")" 2>/dev/null || true
    else
      for sdir in "$combo_dir"/tf/*_r${rid}_l${lid}; do
        [[ -d "$sdir" ]] || continue
        site="$(basename "$sdir")"
        mkdir -p "$CLAIM/$(claim_name "$combo" "$site")" 2>/dev/null || true
      done
    fi
  done < <(pgrep -af 'run_filt.py' || true)
}

claim_one() {
  local combo_dir combo sdir site tfs tfe cname
  for combo_dir in "$ROOT"/*; do
    [[ -d "$combo_dir/met" && -d "$combo_dir/tf" && -d "$combo_dir/filter" ]] || continue
    combo="$(basename "$combo_dir")"
    [[ "$combo" == logs ]] && continue
    for sdir in "$combo_dir"/tf/*_r*_l*; do
      [[ -d "$sdir" ]] || continue
      site="$(basename "$sdir")"
      [[ -f "$combo_dir/filter/$site/wd.out" ]] && continue
      tfs="$sdir/tfs"
      tfe="$sdir/tfe"
      if ! has_rank_csv "$tfs" || ! has_rank_csv "$tfe"; then
        continue
      fi
      cname="$(claim_name "$combo" "$site")"
      if mkdir "$CLAIM/$cname" 2>/dev/null; then
        echo "$combo $site"
        return 0
      fi
    done
  done
  return 1
}

run_claimed() {
  local combo="$1" site="$2"
  local combo_dir="$ROOT/$combo"
  local bench rid loc
  bench="${site%%.*}"
  rid="${site#*_r}"
  rid="${rid%%_l*}"
  loc="${site##*_l}"
  local site_log="$ROOT/logs/filt_${combo}_${site}.log"
  site_log="${site_log//+/_}"
  : > "$site_log"
  {
    echo "[$(date '+%F %T')] ==> FILT $combo $site rid=$rid lid=$loc worker=$$"
  } | flock "$LOG" tee -a "$LOG" >/dev/null
  if ${PYTHON} "$SCRIPT_DIR/run_filt.py" "$combo_dir" \
       --kernel "$bench" --rid "$rid" --lid "$loc" \
       --filt "$SUITE/bin/filt.x" >> "$site_log" 2>&1; then
    {
      echo "[$(date '+%F %T')] OK   FILT $combo $site"
    } | flock "$LOG" tee -a "$LOG" >/dev/null
    return 0
  fi
  {
    echo "[$(date '+%F %T')] FAIL FILT $combo $site (see $site_log)"
  } | flock "$LOG" tee -a "$LOG" >/dev/null
  return 1
}

worker() {
  local wid="$1" pair combo site
  while pair="$(claim_one)"; do
    combo="${pair%% *}"
    site="${pair#* }"
    run_claimed "$combo" "$site" || true
  done
  log "worker $wid idle (no unclaimed sites) ROOT=$ROOT"
}

claim_inflight
log "ROOT=$ROOT filt=$SUITE/bin/filt.x workers=$NWORKERS claim=$CLAIM"

for i in $(seq 1 "$NWORKERS"); do
  worker "$i" &
done
wait

nleft=0
for combo_dir in "$ROOT"/*; do
  [[ -d "$combo_dir/filter" ]] || continue
  for sdir in "$combo_dir"/tf/*_r*_l*; do
    [[ -d "$sdir" ]] || continue
    site="$(basename "$sdir")"
    [[ -f "$combo_dir/filter/$site/wd.out" ]] || nleft=$((nleft + 1))
  done
done
if [[ $nleft -eq 0 ]]; then
  : > "$ROOT/logs/filt_all.done"
  log "filt_all finished. parallel workers done dir=$ROOT"
else
  log "workers exited with $nleft sites still missing wd.out dir=$ROOT"
fi
