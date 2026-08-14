#!/usr/bin/env bash
# =============================================================================
# run_matrix.sh — NPB-MPI TacVar 全矩阵测量驱动（在 c920bn1 上运行）
#
# 矩阵: 8 kernel (bt cg ep ft is lu mg sp) × 6 timer/counter 组合
#    T1 clock_gettime        + none
#    T2 clock_gettime        + papi_read×4
#    T3 cntvct_el0           + none          (TACVAR_NSTP=10)
#    T4 cntvct_el0           + papi_read×4   (TACVAR_NSTP=10)
#    T5 papi_get_real_nsec   + none
#    T6 papi_get_real_nsec   + papi_read×4   (cg 已完成, manifest 预置跳过)
#
# 每个组合: 常规测量 → nspg标定(每timer一次,跨组合复用) → median
# 幂等: manifest.csv 已有记录的组合跳过; 中断后重跑只补缺。
# 节点礼貌: 每组合前检查 load, 超阈值等待; 全程串行; 单组合失败不中断整体。
# 输出: 每组合一行摘要到 stdout + matrix/logs/<tag>.log; 数据在 matrix/data_<tag>/。
# 本脚本不使用任何 git 命令。
#
# 用法: scripts/run_matrix.sh [--only bench1,bench2] [--timer T1,T3] [--class C]
# =============================================================================
set -uo pipefail

# ---------------- 顶部用户参数 ----------------
SUITE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATRIX_DIR="$SUITE/matrix"
LOG_DIR="$MATRIX_DIR/logs"
MANIFEST="$MATRIX_DIR/manifest.csv"
CONF="$SUITE/tacvar.conf"

BENCH_LIST=(bt cg ep ft is lu mg sp)
TIMER_LIST=(T1 T2 T3 T4 T5 T6)
CLASS="${CLASS:-C}"
NP="${NP:-128}"
# BT/SP 都在源码内硬性检查 "Expecting a square number of processes"，
# 128 进程直接 abort（SP 的 WARNING 文本是误导，实测为 abort），用 11×11=121；
# 其余 kernel 128。可在下方按需覆盖。
declare -A BENCH_NP=( [bt]=121 [sp]=121 )
LOAD_MAX="${LOAD_MAX:-2.0}"          # load average 超此值等待（勿打扰节点）。
                                     # 0.5 太严: 自身 121/128 进程测量后的 load 衰减到 0.5
                                     # 需 2-3 分钟, 46 组合空等 ~2h; 2.0 仍能挡住外部大任务。
POLL_INT="${POLL_INT:-60}"           # load 等待轮询秒数

MPI_HOME=/home/hpckey/01-App/openmpi-5.0.8
PAPI_HOME=/home/hpckey/01-App/papi
export PATH="$MPI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_HOME/lib:$PAPI_HOME/lib:${LD_LIBRARY_PATH:-}"

# ---------------- 组合定义 ----------------
# 格式: 标签|TIMER|COUNTER_BACKEND|COUNTER_COUNT|COUNTER_NAMES|NSTP
declare -A COMBO=(
  [T1]="clock_gettime|none|0||0"
  [T2]="clock_gettime|papi_read|4|CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS|0"
  [T3]="cntvct_el0|none|0||10"
  [T4]="cntvct_el0|papi_read|4|CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS|10"
  [T5]="papi_get_real_nsec|none|0||0"
  [T6]="papi_get_real_nsec|papi_read|4|CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS|0"
)

# ---------------- 解析参数 ----------------
ONLY_BENCH=""
ONLY_TIMER=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY_BENCH="$2"; shift 2 ;;
    --timer) ONLY_TIMER="$2"; shift 2 ;;
    --class) CLASS="$2"; shift 2 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done

[[ -n "$ONLY_BENCH" ]] && BENCH_LIST=($(echo "$ONLY_BENCH" | tr ',' ' '))
[[ -n "$ONLY_TIMER" ]] && TIMER_LIST=($(echo "$ONLY_TIMER" | tr ',' ' '))
BENCH_U="${BENCH_LIST[0]^^}"

mkdir -p "$MATRIX_DIR" "$LOG_DIR"

# ---------------- 工具函数 ----------------
log() { echo "[$(date '+%F %T')] $*"; }

set_conf_kv() {   # set_conf_kv KEY VALUE
  local key="$1" val="$2"
  if grep -q "^${key}=" "$CONF"; then
    sed -i "s|^${key}=.*|${key}=${val}|" "$CONF"
  else
    echo "${key}=${val}" >> "$CONF"
  fi
}

wait_for_idle() {  # 节点礼貌: load 超阈值则等待
  while :; do
    local ld
    ld=$(cut -d' ' -f1 /proc/loadavg)
    if awk -v l="$ld" -v m="$LOAD_MAX" 'BEGIN{exit !(l<=m)}'; then
      return 0
    fi
    log "load=${ld} > ${LOAD_MAX}, sleep ${POLL_INT}s (do not disturb node)"
    sleep "$POLL_INT"
  done
}

is_done() {  # is_done BENCH CLASS TAG ; FAILED 行不算完成
  [[ -f "$MANIFEST" ]] || return 1
  grep -q "^${1},${2},${3},.*,SUCCESSFUL,\|^${1},${2},${3},.*,UNVERIFIED," "$MANIFEST"
}

run_mpi() {  # run_mpi EXTRA_ENV... -- BIN [ARGS...]
  local envs=()
  while [[ $# -gt 0 && "$1" != "--" ]]; do
    envs+=("$1"); shift
  done
  [[ $# -gt 0 && "$1" == "--" ]] && shift
  local bin="$1"; shift
  env "${envs[@]}" NPB_TIMER_FLAG=1 \
    mpirun --map-by core --bind-to core -np "$np_this" "$bin" "$@"
}

# ---------------- 主循环 ----------------
# 若 manifest 不存在则初始化表头（并预置已完成组合）
if [[ ! -f "$MANIFEST" ]]; then
  echo "bench,class,tag,timer,counter,data_dir,verification,nspg,median_sites,er_summary,finished_at" > "$MANIFEST"
  # cg.C × T6 已完成于 data_20260813T223213（本轮之前跑过），预置跳过
  echo "cg,C,T6,papi_get_real_nsec,papi_read×4,data_20260813T223213,SUCCESSFUL,0.344825269,12,er见data_20260813T223213,2026-08-13T22:50" >> "$MANIFEST"
  log "manifest initialized: $MANIFEST"
fi

declare -A NSPG_DONE=()   # timer -> nspg.txt 路径（每 timer 只标定一次）

for bench in "${BENCH_LIST[@]}"; do
  bench_u="${bench^^}"
  np_this="${BENCH_NP[$bench]:-$NP}"
  for tag in "${TIMER_LIST[@]}"; do
    IFS='|' read -r timer cbackend ccount cnames nstp <<< "${COMBO[$tag]}"
    [[ -z "$timer" ]] && { echo "unknown combo $tag"; exit 2; }

    if is_done "$bench" "$CLASS" "$tag"; then
      log "skip ${bench}.${CLASS} ${tag} (already in manifest)"
      continue
    fi
    wait_for_idle

    data_dir_abs="$MATRIX_DIR/data_${bench}_${tag}"
    data_dir_rel="matrix/data_${bench}_${tag}"
    run_log="$LOG_DIR/${bench}_${CLASS}_${tag}.log"
    t0=$(date +%s)
    log "==> START ${bench}.${CLASS} ${tag} (${timer} + ${cbackend:-none}, np=${np_this})"
    {
      echo "### run_matrix: ${bench}.${CLASS} ${tag} at $(date)"
    } > "$run_log"

    ok_verification=""
    nspg_val=""
    failed=0

    # ---- 1) 配置 conf: timer/counter, TF=OFF ----
    set_conf_kv TACVAR_TIMER "$timer"
    set_conf_kv TACVAR_COUNTER_BACKEND "$cbackend"
    set_conf_kv TACVAR_COUNTER_COUNT "$ccount"
    set_conf_kv TACVAR_COUNTER_NAMES "$cnames"
    set_conf_kv TACVAR_NSTP "$nstp"
    set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
    set_conf_kv TACVAR_TF_DATA_ROOT ""

    # ---- 2) 构建常规二进制 ----
    if ! (cd "$SUITE" && make tacvar_clean && make "$bench_u" CLASS="$CLASS" PAPI_HOME="$PAPI_HOME") >> "$run_log" 2>&1; then
      log "FAIL ${bench}.${CLASS} ${tag}: build (normal) failed"
      failed=1
    fi

    # ---- 3) nspg 标定：构建（每 timer 一次） ----
    if [[ $failed -eq 0 && -z "${NSPG_DONE[$timer]:-}" ]]; then
      if ! (cd "$SUITE" && make nspg NSPG_MPICC="$(command -v mpicc)" PAPI_HOME="$PAPI_HOME") >> "$run_log" 2>&1; then
        log "FAIL ${bench}.${CLASS} ${tag}: build nspg failed"
        failed=1
      fi
    fi

    # ---- 4) 常规测量 ----
    if [[ $failed -eq 0 ]]; then
      rm -rf "$data_dir_abs"
      if run_mpi TACVAR_DATA_DIR="$data_dir_abs" -- "$SUITE/bin/${bench}.${CLASS}.x" >> "$run_log" 2>&1; then
        if grep -qE 'SUCCESSFUL' "$run_log"; then
          ok_verification=SUCCESSFUL
        else
          ok_verification=UNVERIFIED
          log "WARN ${bench}.${CLASS} ${tag}: no SUCCESSFUL in stdout"
        fi
      else
        log "FAIL ${bench}.${CLASS} ${tag}: kernel run failed"
        failed=1
      fi
    fi

    # ---- 5) nspg 运行（每 timer 一次, 复制到各 data_dir）+ median ----
    if [[ $failed -eq 0 ]]; then
      if [[ -z "${NSPG_DONE[$timer]:-}" ]]; then
        if run_mpi -- "$SUITE/bin/test_nspg.x" "$data_dir_abs/nspg.txt" >> "$run_log" 2>&1; then
          NSPG_DONE[$timer]="$data_dir_abs/nspg.txt"
        else
          log "FAIL ${bench}.${CLASS} ${tag}: nspg run failed"
          failed=1
        fi
      else
        cp "${NSPG_DONE[$timer]}" "$data_dir_abs/nspg.txt"
      fi
      if [[ $failed -eq 0 ]]; then
        nspg_val=$(awk '{print $1}' "$data_dir_abs/nspg.txt" 2>/dev/null | head -1)
        if ! (cd "$SUITE" && python3 scripts/get_median.py "$data_dir_abs") >> "$run_log" 2>&1; then
          log "FAIL ${bench}.${CLASS} ${tag}: get_median failed"
          failed=1
        fi
      fi
    fi

    # ---- 6) 恢复 conf + manifest ----
    set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
    set_conf_kv TACVAR_TF_DATA_ROOT ""
    t1=$(date +%s)
    if [[ $failed -eq 0 ]]; then
      echo "${bench},${CLASS},${tag},${timer},${cbackend},${data_dir_rel},${ok_verification},${nspg_val},?,,$(date '+%F %T')" >> "$MANIFEST"
      log "OK   ${bench}.${CLASS} ${tag}  [$(($t1-$t0))s]  ver=${ok_verification} nspg=${nspg_val}"
    else
      echo "${bench},${CLASS},${tag},${timer},${cbackend},${data_dir_rel},FAILED,${nspg_val},?,,$(date '+%F %T')" >> "$MANIFEST"
      log "FAIL ${bench}.${CLASS} ${tag}  [$(($t1-$t0))s]  see $run_log"
    fi
  done
done

log "matrix finished. manifest: $MANIFEST"
exit 0