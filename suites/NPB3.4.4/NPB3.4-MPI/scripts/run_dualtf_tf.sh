#!/usr/bin/env bash
# =============================================================================
# run_dualtf_tf.sh — 按 median.csv 的 ngauge 编 tfs/tfe，逐站点双采样
#
# 无命令行传参。改 dualtf_common.sh 顶部变量区后执行：
#   bash scripts/run_dualtf_tf.sh
# 写出目录：<combo>/tf/kernel.class_r${rid}_l${loc}/{tfs,tfe}/
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
    med="$combo/met/median.csv"
    if [[ ! -f "$med" ]]; then
      log "WARN skip combo=$(combo_name "$timer" "$reader"): missing $med（先跑 run_dualtf_met.sh）"
      continue
    fi
    mkdir -p "$combo/tf"
    cp -f "$ROOT/nspg.txt" "$combo/tf/nspg.txt"
    cp -f "$med" "$combo/tf/median.csv"
    log "==> TF combo=$(combo_name "$timer" "$reader")"

    for spec in "${KERNELS[@]}"; do
      parse_kernel "$spec"
      np_this="$(np_for_kernel "$spec")"

      while read -r pair; do
        [[ -z "$pair" ]] && continue
        rid="${pair%%:*}"
        loc="${pair##*:}"
        site="${_kc}_r${rid}_l${loc}"
        dst_tfs="$combo/tf/${site}/tfs"
        dst_tfe="$combo/tf/${site}/tfe"
        if has_rank_csv "$dst_tfs" && has_rank_csv "$dst_tfe"; then
          log "skip TF $site (has tfs+tfe CSV)"
          continue
        fi

        site_log="$ROOT/logs/tf_$(combo_name "$timer" "$reader")_${site}.log"
        site_log="${site_log//+/_}"
        : > "$site_log"
        log "==> TF $site timer=$timer reader=$reader np=$np_this"

        apply_combo_conf "$timer" "$reader"
        set_conf_kv TACVAR_TF_SAMPLING_MODE ON
        set_conf_kv TACVAR_TF_DATA_ROOT "$combo/tf"
        set_conf_kv TACVAR_TF_REG_ID "$rid"
        set_conf_kv TACVAR_TF_LOC_ID "$loc"
        set_conf_kv TACVAR_NSPT_FILE "$ROOT/nspt.txt"

        if ! (cd "$SUITE" && make tacvar_clean && make "$_bench_u" CLASS="$_class" PAPI_HOME="$PAPI_HOME") >> "$site_log" 2>&1 \
           || [[ ! -x "$SUITE/bin/${_kc}_tfs.x" ]] || [[ ! -x "$SUITE/bin/${_kc}_tfe.x" ]]; then
          log "FAIL TF $site: build (see $site_log)"
          continue
        fi

        rm -rf "$combo/tf/${_kc}_tfs" "$combo/tf/${_kc}_tfe"
        tfs_ok=0
        tfe_ok=0
        if run_mpi "$combo/tf" "$np_this" "$SUITE/bin/${_kc}_tfs.x" > "$site_log.tmp" 2>&1 \
           && grep -q SUCCESSFUL "$site_log.tmp"; then
          tfs_ok=1
          cat "$site_log.tmp" >> "$site_log"
        else
          cat "$site_log.tmp" >> "$site_log" 2>/dev/null || true
          log "FAIL TF $site: tfs not SUCCESSFUL"
        fi
        if [[ $tfs_ok -eq 1 ]] \
           && run_mpi "$combo/tf" "$np_this" "$SUITE/bin/${_kc}_tfe.x" > "$site_log.tmp" 2>&1 \
           && grep -q SUCCESSFUL "$site_log.tmp"; then
          tfe_ok=1
          cat "$site_log.tmp" >> "$site_log"
        else
          [[ $tfs_ok -eq 1 ]] && {
            cat "$site_log.tmp" >> "$site_log" 2>/dev/null || true
            log "FAIL TF $site: tfe not SUCCESSFUL"
          }
        fi
        rm -f "$site_log.tmp"

        if [[ $tfs_ok -eq 1 && -d "$combo/tf/${_kc}_tfs" ]]; then
          rm -rf "$dst_tfs"
          mkdir -p "$(dirname "$dst_tfs")"
          mv "$combo/tf/${_kc}_tfs" "$dst_tfs"
        fi
        if [[ $tfe_ok -eq 1 && -d "$combo/tf/${_kc}_tfe" ]]; then
          rm -rf "$dst_tfe"
          mkdir -p "$(dirname "$dst_tfe")"
          mv "$combo/tf/${_kc}_tfe" "$dst_tfe"
        fi
        if [[ $tfs_ok -eq 1 && $tfe_ok -eq 1 ]]; then
          log "OK   TF $site"
        fi
      done < <(sites_for_kernel "$_bench")
    done
  done
done

set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
set_conf_kv TACVAR_TF_DATA_ROOT ""
log "tf finished. ROOT=$ROOT"
exit 0
