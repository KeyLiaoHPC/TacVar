#!/usr/bin/env bash
# =============================================================================
# run_dualtf_filt_all.sh — cicx6330n1 only: filt every complete combo × site
#
# Discovers combo directory names under ROOT (cntvct_el0+… or rdtscp+…).
# Calls run_filt.py once per (combo, kernel, rid, lid). Never mixes sites.
# Requires check_dualtf_complete.sh to pass first.
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dualtf_common.sh
source "$SCRIPT_DIR/dualtf_common.sh"

resolve_root || exit 2

if ! bash "$SCRIPT_DIR/check_dualtf_complete.sh" "$ROOT"; then
  log "refuse filt: tree not complete ($ROOT)"
  exit 2
fi

if [[ ! -x "$SUITE/bin/filt.x" ]]; then
  log "missing $SUITE/bin/filt.x; run make filt on cicx6330n1"
  exit 2
fi

log "ROOT=$ROOT filt=$SUITE/bin/filt.x"
nsite=0
nfail=0
nskip=0

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

for combo_dir in "$ROOT"/*; do
  [[ -d "$combo_dir/met" && -d "$combo_dir/tf" && -d "$combo_dir/filter" ]] || continue
  combo="$(basename "$combo_dir")"
  [[ "$combo" == logs ]] && continue
  log "==> FILT combo=$combo"

  for spec in "${KERNELS[@]}"; do
    parse_kernel "$spec"
    while read -r pair; do
      [[ -z "$pair" ]] && continue
      rid="${pair%%:*}"
      loc="${pair##*:}"
      site="${_kc}_r${rid}_l${loc}"
      tfs="$combo_dir/tf/${site}/tfs"
      tfe="$combo_dir/tf/${site}/tfe"
      if ! has_rank_csv "$tfs" || ! has_rank_csv "$tfe"; then
        log "skip FILT $combo $site (incomplete tfs/tfe)"
        nskip=$((nskip + 1))
        continue
      fi
      site_log="$ROOT/logs/filt_${combo}_${site}.log"
      site_log="${site_log//+/_}"
      : > "$site_log"
      log "==> FILT $combo $site rid=$rid lid=$loc"
      if ${PYTHON} "$SCRIPT_DIR/run_filt.py" "$combo_dir" \
           --kernel "$_bench" --rid "$rid" --lid "$loc" \
           --filt "$SUITE/bin/filt.x" >> "$site_log" 2>&1; then
        nsite=$((nsite + 1))
        log "OK   FILT $combo $site"
      else
        nfail=$((nfail + 1))
        log "FAIL FILT $combo $site (see $site_log)"
      fi
    done < <(sites_for_kernel "$_bench")
  done
done

log "filt_all finished. ok=$nsite fail=$nfail skip=$nskip dir=$ROOT"
[[ $nfail -eq 0 && $nsite -gt 0 ]]
exit $?
