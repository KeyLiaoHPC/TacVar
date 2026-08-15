#!/usr/bin/env bash
# =============================================================================
# filt_watch.sh — 本机（x86 工作台）常驻 FilT watcher
#   扫描 matrix/ext_<stamp>/T*/manifest.csv 中 status=MEASURED 且未过滤的站点，
#   用本机编译的 filt.x86.x 直接过滤 NFS 上的数据（避免 aarch64 filt.x 跨架构）
#   （Exec format error：filt.x 在 c920bn1 构建，不能在本机 x86 跑，见技能 pitfall）。
# 用法：
#   bash scripts/filt_watch.sh            # 常驻循环（默认 60s 一轮）
#   bash scripts/filt_watch.sh --once     # 单轮后退出（冒烟/手动）
#   环境变量: FILT_WATCH_INTERVAL（默认 60）、FILT_BIN（默认 bin/filt.x86.x）
# 幂等：已滤站点记录在 <root>/filter_done.tsv，失败重试 3 次后记 filter_fail.tsv。
# =============================================================================
set -uo pipefail
SUITE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FILT_BIN="${FILT_BIN:-$SUITE/bin/filt.x86.x}"
INTERVAL="${FILT_WATCH_INTERVAL:-60}"
ONCE=0
[[ "${1:-}" == "--once" ]] && ONCE=1

# 单实例锁（同主机同时只跑一个 watcher）
LOCK="/tmp/tacvar-filt-watch-${USER}.lock"
exec 9>"$LOCK" || exit 4
flock -n 9 || { echo "[$(date '+%F %T')] another watcher running, exit"; exit 0; }

[[ -x "$FILT_BIN" ]] || { echo "filt.x86.x not found/not executable: $FILT_BIN"; exit 2; }

process_one() {   # process_one BENCH SITE DATA_DIR ROOT
  local bench="$1" site="$2" data_dir="$3" root="$4"
  local rid loc
  rid="${site#r}"; rid="${rid%_*}"; loc="${site#*_}"; loc="${loc#l}"
  local done_file="$root/filter_done.tsv" fail_file="$root/filter_fail.tsv"
  [[ -f "$done_file" ]] && grep -q "^${bench},${site}," "$done_file" 2>/dev/null && return 0
  local err_log="$root/filt_err_${bench}_${site}.log"
  if python3 "$SUITE/scripts/run_filt.py" "$data_dir" --rid "$rid" --lid "$loc" \
      --filt "$FILT_BIN" >"$err_log" 2>&1; then
    local er="" ep=""
    [[ -f "$data_dir/${bench}.C_filt/r${rid}_l${loc}/er.out" ]] && \
      er=$(tr -d '\n' < "$data_dir/${bench}.C_filt/r${rid}_l${loc}/er.out")
    [[ -f "$data_dir/${bench}.C_filt/r${rid}_l${loc}/ep.out" ]] && \
      ep=$(tr -d '\n' < "$data_dir/${bench}.C_filt/r${rid}_l${loc}/ep.out")
    echo -e "$bench\t$site\t${er:-}\t${ep:-}\t$(date '+%F %T')" >> "$done_file"
    echo "[$(date '+%F %T')] filtered $bench $site er=${er:-} ep=${ep:-}"
  else
    local n
    n=$(grep -c "^${bench},${site}," "$fail_file" 2>/dev/null || echo 0)
    echo -e "$bench\t$site\t$(date '+%F %T')" >> "$fail_file"
    if [[ $n -ge 2 ]]; then
      echo "[$(date '+%F %T')] FILT-FAIL(3x) $bench $site see $err_log"
    else
      echo "[$(date '+%F %T')] filt retry $bench $site (attempt $((n+1)))"
    fi
  fi
}

while :; do
  touched=0
  for mf in "$SUITE"/matrix/ext_*/T*/manifest.csv; do
    [[ -f "$mf" ]] || continue
    root="$(dirname "$mf")"
    # status=MEASURED 且不在 done 账本的行 -> 过滤
    awk -F, '$3=="MEASURED"' "$mf" | while IFS= read -r ln; do
      bench="${ln%%,*}"; rest="${ln#*,}"; site="${rest%%,*}"
      [[ "$site" == "site" ]] && continue
      # manifest 列: bench,site,status,...,data_dir(第11列),finished_at(第12列)
      data_dir="$(echo "$ln" | cut -d, -f11)"
      [[ -d "$data_dir" ]] || { echo "[$(date '+%F %T')] skip missing data_dir: $data_dir"; continue; }
      process_one "$bench" "$site" "$data_dir" "$root"
    done
    touched=1
  done
  [[ $ONCE -eq 1 ]] && break
  sleep "$INTERVAL"
done
exit 0