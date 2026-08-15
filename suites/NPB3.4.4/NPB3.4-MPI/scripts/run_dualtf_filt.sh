#!/usr/bin/env bash
# =============================================================================
# run_dualtf_filt.sh — 按 KERNELS × SITES_* 对指定 combo 做 FilT
#
# 无命令行传参。改本文件顶部 COMBO，以及 dualtf_common.sh 的 KERNELS / SITES_*：
#   bash scripts/run_dualtf_filt.sh
# 只滤一站时：把 KERNELS 和对应 SITES_* 收成单项（可在 source 之后覆盖）。
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ===== 本脚本变量 =====
COMBO="native+none"    # 含 met/tf/filter/ 的 combo 目录名（相对 ROOT）
# ========================

# shellcheck source=dualtf_common.sh
source "$SCRIPT_DIR/dualtf_common.sh"
# 只滤子集时在此覆盖，例如：
# KERNELS=(bt.C)
# SITES_bt=(7:1)

resolve_root || exit 2
combo_dir="$ROOT/$COMBO"
if [[ ! -d "$combo_dir/met" || ! -d "$combo_dir/tf" || ! -d "$combo_dir/filter" ]]; then
  log "COMBO 目录须已含 met/ tf/ filter/: $combo_dir"
  exit 2
fi
if [[ ! -x "$SUITE/bin/filt.x" ]]; then
  log "missing $SUITE/bin/filt.x；先跑 run_dualtf_init.sh"
  exit 2
fi

log "ROOT=$ROOT COMBO=$COMBO"
nsite=0
nfail=0

for spec in "${KERNELS[@]}"; do
  parse_kernel "$spec"
  while read -r pair; do
    [[ -z "$pair" ]] && continue
    rid="${pair%%:*}"
    loc="${pair##*:}"
    site="${_kc}_r${rid}_l${loc}"
    site_log="$ROOT/logs/filt_${COMBO}_${site}.log"
    site_log="${site_log//+/_}"
    : > "$site_log"
    log "==> FILT $COMBO $site"
    if ${PYTHON} "$SCRIPT_DIR/run_filt.py" "$combo_dir" \
         --rid "$rid" --lid "$loc" --filt "$SUITE/bin/filt.x" >> "$site_log" 2>&1; then
      nsite=$((nsite + 1))
      log "OK   FILT $site"
    else
      nfail=$((nfail + 1))
      log "FAIL FILT $site (see $site_log)"
    fi
  done < <(sites_for_kernel "$_bench")
done

log "filt finished. ok=$nsite fail=$nfail dir=$combo_dir"
[[ $nfail -eq 0 && $nsite -gt 0 ]]
exit $?
