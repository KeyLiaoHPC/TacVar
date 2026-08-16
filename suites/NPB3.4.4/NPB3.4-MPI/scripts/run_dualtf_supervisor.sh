#!/usr/bin/env bash
# =============================================================================
# run_dualtf_supervisor.sh — cicx6330n1: 5-min NFS report + filt when complete
#
#   bash scripts/run_dualtf_supervisor.sh
#
# Does not SSH. Starts run_dualtf_filt_all.sh once per finished tree.
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WS="/mnt/keylabmain/nfs/hpckey/03-Project/TacVar_NPC"
INTERVAL=300
HOSTS=(c920bn1 cgnr6760pn1 camd9755n2)

log() { echo "[$(date '+%F %T')] $*"; }

while true; do
  bash "$SCRIPT_DIR/watch_dualtf.sh"
  for h in "${HOSTS[@]}"; do
    tree="$WS/TacVar_NPC_${h}/suites/NPB3.4.4/NPB3.4-MPI"
    root="$tree/matrix/${h}_dualtf_20260815"
    mark="$root/logs/filt_all.started"
    donef="$root/logs/filt_all.done"
    [[ -d "$root" ]] || continue
    [[ -f "$donef" ]] && continue
    [[ -f "$mark" ]] && continue
    if bash "$tree/scripts/check_dualtf_complete.sh" "$root" >/dev/null 2>&1; then
      log "complete $h — build filt.x and start filt_all"
      mkdir -p "$root/logs"
      if ! (cd "$tree" && make filt) >> "$root/logs/filt_all.log" 2>&1; then
        log "FAIL make filt for $h"
        continue
      fi
      : > "$mark"
      nohup bash "$tree/scripts/run_dualtf_filt_all.sh" \
        >> "$root/logs/filt_all.log" 2>&1 < /dev/null &
      echo $! > "$root/logs/filt_all.pid"
      log "filt_all pid=$(cat "$root/logs/filt_all.pid") host=$h"
    fi
    if [[ -f "$root/logs/filt_all.pid" ]]; then
      pid="$(cat "$root/logs/filt_all.pid")"
      if [[ -n "$pid" ]] && ! kill -0 "$pid" 2>/dev/null; then
        if grep -q 'filt_all finished' "$root/logs/filt_all.log" 2>/dev/null; then
          : > "$donef"
          log "filt_all done host=$h"
        fi
      fi
    fi
  done
  sleep "$INTERVAL"
done
