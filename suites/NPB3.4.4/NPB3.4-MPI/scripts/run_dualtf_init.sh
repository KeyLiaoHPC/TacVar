#!/usr/bin/env bash
# =============================================================================
# run_dualtf_init.sh — 建 DualTF 目录树，一次性标定 nspt / nspg，编译 filt.x
#
# 无命令行传参。改 dualtf_common.sh 顶部变量区后执行：
#   bash scripts/run_dualtf_init.sh
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dualtf_common.sh
source "$SCRIPT_DIR/dualtf_common.sh"

ensure_root_new
log "ROOT=$ROOT"

# ---- 目录骨架：每个 combo 的 met/ + 每个站点的 tf/ filter/ ----
for timer in "${TIMERS[@]}"; do
  for reader in "${EVENT_READERS[@]}"; do
    combo="$ROOT/$(combo_name "$timer" "$reader")"
    mkdir -p "$combo/met" "$combo/tf" "$combo/filter"
    for spec in "${KERNELS[@]}"; do
      parse_kernel "$spec"
      mkdir -p "$combo/met/${_kc}"
      while read -r pair; do
        [[ -z "$pair" ]] && continue
        rid="${pair%%:*}"
        loc="${pair##*:}"
        mkdir -p "$combo/tf/${_kc}_r${rid}_l${loc}/tfs"
        mkdir -p "$combo/tf/${_kc}_r${rid}_l${loc}/tfe"
        mkdir -p "$combo/filter/${_kc}_r${rid}_l${loc}"
      done < <(sites_for_kernel "$_bench")
    done
    log "tree $combo"
  done
done

init_log="$ROOT/logs/init.log"
: > "$init_log"

# ---- nspt：probe_counter_freq --nspt-out（不是旧的 probe_clock_rate.c）----
if [[ -f "$ROOT/nspt.txt" && "$FORCE_CALIB" -eq 0 ]]; then
  log "skip nspt (exists $ROOT/nspt.txt)"
else
  log "nspt -> $ROOT/nspt.txt"
  if ! (cd "$SUITE" && make probe_counter_freq) >> "$init_log" 2>&1; then
    log "FAIL: make probe_counter_freq (see $init_log)"
    exit 1
  fi
  if ! "$SUITE/common/probe_counter_freq.x" --nspt-out "$ROOT/nspt.txt" >> "$init_log" 2>&1; then
    log "FAIL: probe_counter_freq --nspt-out"
    exit 1
  fi
fi

# ---- nspg：节点 tick 计时器标定一次，全矩阵复用 ----
if [[ -f "$ROOT/nspg.txt" && "$FORCE_CALIB" -eq 0 ]]; then
  log "skip nspg (exists $ROOT/nspg.txt)"
else
  log "nspg timer=$TICK_TIMER -> $ROOT/nspg.txt"
  apply_combo_conf "$TICK_TIMER" none
  set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
  set_conf_kv TACVAR_TF_DATA_ROOT ""
  set_conf_kv TACVAR_NSPT_FILE "$ROOT/nspt.txt"
  if ! (cd "$SUITE" && make tacvar_clean && make nspg NSPG_MPICC="$(command -v mpicc)" PAPI_HOME="$PAPI_HOME") >> "$init_log" 2>&1; then
    log "FAIL: make nspg (see $init_log)"
    exit 1
  fi
  np_nspg="$NP"
  [[ "$np_nspg" -eq 0 ]] && np_nspg="${NODE_NP[$NODE]}"
  if ! NPB_TIMER_FLAG=1 timeout "$RUN_TIMEOUT" \
       mpirun "${MPIRUN_ARGS[@]}" "$np_nspg" \
       "$SUITE/bin/test_nspg.x" "$ROOT/nspg.txt" >> "$init_log" 2>&1; then
    log "FAIL: test_nspg.x (see $init_log)"
    exit 1
  fi
fi

# ---- filt.x ----
if ! (cd "$SUITE" && make filt) >> "$init_log" 2>&1; then
  log "FAIL: make filt (see $init_log)"
  exit 1
fi

log "init done. ROOT=$ROOT"
log "  nspt=$(tr -d '\n' < "$ROOT/nspt.txt" | awk '{print $1}')"
log "  nspg=$(tr -d '\n' < "$ROOT/nspg.txt" | awk '{print $1}')"
exit 0
