#!/usr/bin/env bash
# =============================================================================
# run_dualtf_met.sh — 遍历 timer × event_reader × kernel，常规计时（不采样）
#
# 无命令行传参。改 dualtf_common.sh 顶部变量区后执行：
#   bash scripts/run_dualtf_met.sh
# 每个 kernel 跑完后 sleep，并对该 combo 的 met/ 重算 median.csv / met_stat.csv。
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=dualtf_common.sh
source "$SCRIPT_DIR/dualtf_common.sh"

resolve_root || exit 2
log "ROOT=$ROOT"

if [[ ! -f "$ROOT/nspg.txt" ]]; then
  log "missing $ROOT/nspg.txt；先跑 run_dualtf_init.sh"
  exit 2
fi

for timer in "${TIMERS[@]}"; do
  for reader in "${EVENT_READERS[@]}"; do
    combo="$ROOT/$(combo_name "$timer" "$reader")"
    mkdir -p "$combo/met" "$ROOT/logs"
    cp -f "$ROOT/nspg.txt" "$combo/met/nspg.txt"
    log "==> MET combo=$(combo_name "$timer" "$reader")"

    for spec in "${KERNELS[@]}"; do
      parse_kernel "$spec"
      np_this="$(np_for_kernel "$spec")"
      kdir="$combo/met/${_kc}"
      if has_rank_csv "$kdir"; then
        log "skip MET ${_kc} (has CSV in $kdir)"
        continue
      fi

      bench_log="$ROOT/logs/met_$(combo_name "$timer" "$reader")_${_kc}.log"
      bench_log="${bench_log//+/_}"
      : > "$bench_log"
      log "==> MET ${_kc} timer=$timer reader=$reader np=$np_this"

      apply_combo_conf "$timer" "$reader"
      set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
      set_conf_kv TACVAR_TF_DATA_ROOT ""
      set_conf_kv TACVAR_NSPT_FILE "$ROOT/nspt.txt"

      if ! (cd "$SUITE" && make tacvar_clean && make "$_bench_u" CLASS="$_class" PAPI_HOME="$PAPI_HOME") >> "$bench_log" 2>&1; then
        log "FAIL MET ${_kc}: build (see $bench_log)"
        continue
      fi

      rm -rf "$kdir"
      if run_mpi "$combo/met" "$np_this" "$SUITE/bin/${_kc}.x" > "$bench_log.tmp" 2>&1 \
         && grep -q SUCCESSFUL "$bench_log.tmp"; then
        cat "$bench_log.tmp" >> "$bench_log"
        rm -f "$bench_log.tmp"
        if ! ${PYTHON} "$SCRIPT_DIR/get_met_stat.py" "$combo/met" --nspg "$ROOT/nspg.txt" >> "$bench_log" 2>&1; then
          log "FAIL MET ${_kc}: get_met_stat.py"
          continue
        fi
        log "OK   MET ${_kc}"
      else
        cat "$bench_log.tmp" >> "$bench_log" 2>/dev/null || true
        rm -f "$bench_log.tmp"
        log "FAIL MET ${_kc}: run not SUCCESSFUL (see $bench_log)"
      fi
    done
  done
done

set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
set_conf_kv TACVAR_TF_DATA_ROOT ""
log "met finished. ROOT=$ROOT"
exit 0
