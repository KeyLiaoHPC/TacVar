#!/usr/bin/env bash
# =============================================================================
# dualtf_common.sh — DualTF 四脚本共享配置与函数（无命令行传参）
#
# 改范围只编辑下面「可编辑变量区」。各 runner：
#   bash scripts/run_dualtf_init.sh
#   bash scripts/run_dualtf_met.sh
#   bash scripts/run_dualtf_tf.sh
#   bash scripts/run_dualtf_filt.sh
# 本文件由上述脚本 source，不要直接执行。
# =============================================================================

# ===== 可编辑变量区 =====
NODE="$(hostname -s)"          # camd9755n2 的 hostname 不是节点名，须改成 camd9755n2
STAMP="$(date +%Y%m%dT%H%M%S)" # init 在 ROOT 为空时用此戳建目录
ROOT=""                        # 空：init 新建；met/tf/filt 读 matrix/latest_dualtf.txt
NP=0                           # 0 = 用节点表默认（128/128/256）
SLEEP_AFTER_KERNEL=5
RUN_TIMEOUT=1500               # 单次 mpirun 超时（秒）
FORCE_CALIB=0                  # 1 = init 重算 nspt/nspg
PYTHON="${PYTHON:-python3}"

# x86 节点会把 cntvct_el0 自动换成 rdtscp（见 normalize_timers）
TIMERS=(native cntvct_el0 papi_get_real_nsec)
EVENT_READERS=(none papi_read)
KERNELS=(bt.C cg.C ep.C ft.C lu.C mg.C sp.C)

# rid/lid：kernel 名（小写）-> "rid:lid" 列表（2026-08-16 三节点 45 站）
SITES_bt=(7:1 7:2 8:1 9:1 9:2 10:1 10:2 11:1 11:2)
SITES_cg=(3:2 3:3 3:4 3:5)
SITES_ep=(2:1 3:1)
SITES_ft=(6:1 6:2 7:1 7:2 7:3 7:4)
SITES_lu=(3:1 4:1 7:1 7:2 8:1 8:2 9:2)
SITES_mg=(4:1 5:1 6:1 8:3 8:4 8:5 8:6 8:7)
SITES_sp=(6:1 6:2 7:1 8:2 8:5 9:2 9:5 10:2 10:5)

# MPI：${np} 由 np_for_kernel 填入。实际：mpirun "${MPIRUN_ARGS[@]}" "$np" "$bin"
MPIRUN_ARGS=(--map-by core --bind-to core -np)
# ========================

# ---- 以下一般不用改：节点表与函数 ----
SUITE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONF="$SUITE/tacvar.conf"
LATEST_FILE="$SUITE/matrix/latest_dualtf.txt"
NSTP=10

declare -A NODE_MPI_HOME=(
  [c920bn1]=/home/hpckey/01-App/openmpi-5.0.8
  [cgnr6760pn1]=/home/hpckey/01-App/openmpi-5.0.9
  [camd9755n2]=/home/hpckey/01-App/openmpi-5.0.8
)
declare -A NODE_MPI_LIB=(
  [c920bn1]=/home/hpckey/01-App/openmpi-5.0.8/lib
  [cgnr6760pn1]=/home/hpckey/01-App/openmpi-5.0.9/lib
  [camd9755n2]=/home/hpckey/01-App/openmpi-5.0.8/lib
)
declare -A NODE_PAPI_HOME=(
  [c920bn1]=/home/hpckey/01-App/papi
  [cgnr6760pn1]=/home/hpckey/01-App/papi
  [camd9755n2]=/home/hpckey/01-App/papi
)
declare -A NODE_TICK=(
  [c920bn1]=cntvct_el0 [cgnr6760pn1]=rdtscp [camd9755n2]=rdtscp
)
declare -A NODE_NP=(
  [c920bn1]=128 [cgnr6760pn1]=128 [camd9755n2]=256
)
declare -A NODE_PAPI_EVENTS=(
  [c920bn1]="CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS"
  [cgnr6760pn1]="perf::PERF_COUNT_HW_CPU_CYCLES,perf::PERF_COUNT_HW_INSTRUCTIONS,perf::L1-DCACHE-LOADS,perf::PERF_COUNT_HW_BRANCH_INSTRUCTIONS"
  [camd9755n2]="perf::PERF_COUNT_HW_CPU_CYCLES,perf::PERF_COUNT_HW_INSTRUCTIONS,perf::L1-DCACHE-LOADS,perf::PERF_COUNT_HW_BRANCH_INSTRUCTIONS"
)

MPI_HOME="${NODE_MPI_HOME[$NODE]:?未知节点: $NODE (允许 c920bn1/cgnr6760pn1/camd9755n2)}"
MPI_LIB="${NODE_MPI_LIB[$NODE]:-$MPI_HOME/lib}"
PAPI_HOME="${NODE_PAPI_HOME[$NODE]}"
TICK_TIMER="${NODE_TICK[$NODE]}"
PAPI_EVENTS="${NODE_PAPI_EVENTS[$NODE]}"
export PATH="$MPI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_LIB:$PAPI_HOME/lib:${LD_LIBRARY_PATH:-}"
export CC="$(command -v mpicc)"

log() { echo "[$(date '+%F %T')] $*"; }

# 与 NP 距离最小的平方数；并列时取不超过 NP 的那个。
nearest_square_np() {
  local n="$1"
  awk -v n="$n" 'BEGIN{
    if (n < 1) { print 1; exit }
    f = int(sqrt(n));
    if (f < 1) f = 1;
    c = f + 1;
    fs = f * f;
    cs = c * c;
    df = n - fs;
    dc = cs - n;
    if (dc < df) print cs;
    else print fs;
  }'
}

# bt/sp 取最近平方进程数；其它 kernel 用 NP（0 则节点默认）。
np_for_kernel() {
  local spec="$1"
  local bench="${spec%%.*}"
  bench="${bench,,}"
  local base="$NP"
  if [[ "$base" -eq 0 ]]; then
    base="${NODE_NP[$NODE]}"
  fi
  if [[ "$bench" == "bt" || "$bench" == "sp" ]]; then
    nearest_square_np "$base"
  else
    echo "$base"
  fi
}

# 把 TIMERS 里的架构 tick 名对齐到本节点（cntvct_el0 <-> rdtscp）。
normalize_timers() {
  local i t
  for i in "${!TIMERS[@]}"; do
    t="${TIMERS[$i]}"
    if [[ "$t" == "cntvct_el0" && "$TICK_TIMER" == "rdtscp" ]]; then
      TIMERS[$i]=rdtscp
    elif [[ "$t" == "rdtscp" && "$TICK_TIMER" == "cntvct_el0" ]]; then
      TIMERS[$i]=cntvct_el0
    fi
  done
}

# 解析 kernel.class -> _bench _class _bench_u _kc
parse_kernel() {
  local spec="$1"
  _bench="${spec%%.*}"
  _class="${spec#*.}"
  if [[ "$_bench" == "$_class" ]]; then
    _class="C"
  fi
  _bench="${_bench,,}"
  _class="${_class^^}"
  _bench_u="${_bench^^}"
  _kc="${_bench}.${_class}"
}

# 打印该 kernel 的 "rid:lid" 行。
sites_for_kernel() {
  local bench="$1"
  local var="SITES_${bench}"
  declare -p "$var" &>/dev/null || return 0
  local -n _sites="$var"
  [[ ${#_sites[@]} -eq 0 ]] && return 0
  printf '%s\n' "${_sites[@]}"
}

combo_name() { echo "$1+$2"; }

set_conf_kv() {
  local key="$1" val="$2"
  if grep -q "^${key}=" "$CONF"; then
    sed -i "s|^${key}=.*|${key}=${val}|" "$CONF"
  else
    echo "${key}=${val}" >> "$CONF"
  fi
}

# 按 timer × event_reader 写 tacvar.conf 的计时/计数段。
apply_combo_conf() {
  local timer="$1" reader="$2"
  local count=0 names=""
  if [[ "$reader" != "none" && -n "$PAPI_EVENTS" ]]; then
    names="$PAPI_EVENTS"
    count=$(awk -F, '{print NF}' <<< "$names")
  fi
  set_conf_kv TACVAR_TIMER "$timer"
  set_conf_kv TACVAR_NSTP "$NSTP"
  set_conf_kv TACVAR_COUNTER_BACKEND "$reader"
  set_conf_kv TACVAR_COUNTER_COUNT "$count"
  set_conf_kv TACVAR_COUNTER_NAMES "$names"
  if [[ -n "${ROOT:-}" && -f "$ROOT/nspt.txt" ]]; then
    set_conf_kv TACVAR_NSPT_FILE "$ROOT/nspt.txt"
  fi
}

# init：ROOT 空则新建。其它阶段：空则读 latest。
ensure_root_new() {
  mkdir -p "$SUITE/matrix"
  if [[ -z "${ROOT:-}" ]]; then
    ROOT="$SUITE/matrix/${NODE}_dualtf_${STAMP}"
  fi
  mkdir -p "$ROOT/logs"
  echo "$ROOT" > "$LATEST_FILE"
}

resolve_root() {
  if [[ -n "${ROOT:-}" ]]; then
    if [[ ! -d "$ROOT" ]]; then
      log "ROOT 不是目录: $ROOT"
      return 1
    fi
    return 0
  fi
  if [[ -f "$LATEST_FILE" ]]; then
    ROOT="$(tr -d '\n' < "$LATEST_FILE")"
  fi
  if [[ -z "${ROOT:-}" || ! -d "$ROOT" ]]; then
    log "ROOT 为空且没有 $LATEST_FILE；先跑 run_dualtf_init.sh"
    return 1
  fi
}

# run_mpi DATA_DIR NP BIN [args...]；结束后 sleep SLEEP_AFTER_KERNEL。
run_mpi() {
  local data_dir="$1" np="$2" bin="$3"
  shift 3
  NPB_TIMER_FLAG=1 TACVAR_DATA_DIR="$data_dir" \
    timeout "$RUN_TIMEOUT" mpirun "${MPIRUN_ARGS[@]}" "$np" "$bin" "$@" < /dev/null
  local rc=$?
  sleep "$SLEEP_AFTER_KERNEL"
  return $rc
}

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

normalize_timers
