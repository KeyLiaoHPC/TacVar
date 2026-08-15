#!/usr/bin/env bash
# =============================================================================
# run_matrix_insitu.sh — NPB-MPI TacVar 原位双采样矩阵驱动（在 c920bn1 上运行）
#
# 背景：与 2026-08-13 旧单点流程（run_tf_filt.py）不同，新流程在每个站点
#   (region_id, loc_id) 单独构建并运行两次：
#     tfs = 前采样（timer_start 之后、kernel 之前插入 subtraction gauge）
#     tfe = 后采样（kernel 之后、timer_stop 之前插入 subtraction gauge）
#   tf = (tfs - ngauge*nspg) + (tfe' - ngauge*nspg)，tfe 打乱后配对（run_filt.py,
#   seed=1 可复现）；met 取自常规计时跑（Kernel.C/）；migrated=1 行在过滤前丢弃。
# 站点资格（双重门，按本轮实测 met 数据复核，不只看归档）：
#   1) 每个 rank 测点数 >= MIN_SAMPLES_PER_RANK；
#   2) 站点中位数 p50 < MEDIAN_NS_MAX（1 ms）。
# 幂等：manifest.csv（SUCCESSFUL 行跳过）；FAILED 行续跑重试；全程串行；
#   load 门控（不打扰节点）；本脚本不使用任何 git 命令。
#
# 用法（在 c920bn1 上执行，NFS 与工作台共享同一路径）：
#   bash scripts/run_matrix_insitu.sh                          # 全矩阵（77 站）
#   bash scripts/run_matrix_insitu.sh --only bt                # 单 kernel
#   bash scripts/run_matrix_insitu.sh --only bt --site r9_l2   # 单站冒烟
#   bash scripts/run_matrix_insitu.sh --combo T1,T2,T4,T5,T6 --sites good \
#        --no-wait --no-filt                                   # 扩展：计时器×事件
#   --combo Tn,Tm: 组合循环（独立根 matrix/ext_<stamp>/T<t>/）；缺省=单组合兼容旧行为
#   --sites good: 只用 T3 实测 er<0.05 的 49 个良质站点；缺省=77 站全集
#   --no-wait: 关闭 load 门控（独占串行）；--no-filt: 采样后不跑 FilT（本机 watcher 处理）
# =============================================================================
set -uo pipefail

# ---------------- 顶部用户参数 ----------------
SUITE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_LIST=(bt cg ep ft is lu mg sp)          # 矩阵 kernel
CLASS="C"
NP=128                                        # 默认进程数
declare -A BENCH_NP=( [bt]=121 [sp]=121 )     # BT/SP 源码硬性要求平方进程数
LOAD_MAX="${LOAD_MAX:-2.0}"                   # load 超过此值等待（勿打扰节点）
POLL_INT="${POLL_INT:-60}"                    # load 等待轮询秒数
MIN_SAMPLES_PER_RANK=10                       # 单 rank 最少测点数
MEDIAN_NS_MAX=1000000                         # 站点中位数上限：1 ms（ns）
RUN_TIMEOUT=1500                              # 单次 mpirun 超时上限（秒）
STAMP="${STAMP:-$(date +%Y%m%dT%H%M%S)}"      # 本轮矩阵戳（可覆盖）

MPI_HOME=/home/hpckey/01-App/openmpi-5.0.8
PAPI_HOME=/home/hpckey/01-App/papi
export PATH="$MPI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_HOME/lib:$PAPI_HOME/lib:${LD_LIBRARY_PATH:-}"

# 计时组合默认值（单组合模式；扩展模式由 COMBO_TAGS 覆盖）
TIMER=cntvct_el0
COUNTER_BACKEND=none
COUNTER_COUNT=0
COUNTER_NAMES=""
NSTP=10

# 计时×事件组合表（与 2026-08-13 旧矩阵 run_matrix.sh 同源）
# 格式: 标签|TIMER|COUNTER_BACKEND|COUNTER_COUNT|COUNTER_NAMES|NSTP
# NSTP 恒为 10：aarch64 上 tfs/tfe 的 gauge tick 时钟总是 cntvct_el0（10ns/tick，
# tacvar_tf.c 的 TICK_TO_NS 用 NSTP 换算；非 tick 计时器经 per-CPU offset 对齐），
# 与 TACVAR_TIMER 无关。NSTP<10 会让 tikc 侧 elapsed 偏小（冒烟实测 er 爆炸）。
declare -A COMBO_TAGS=(
  [T1]="clock_gettime|none|0||10"
  [T2]="clock_gettime|papi_read|4|CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS|10"
  [T3]="cntvct_el0|none|0||10"
  [T4]="cntvct_el0|papi_read|4|CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS|10"
  [T5]="papi_get_real_nsec|none|0||10"
  [T6]="papi_get_real_nsec|papi_read|4|CPU_CYCLES,INST_RETIRED,PAPI_STL_ICY,PAPI_LST_INS|10"
)

# ---------------- 站点表（严格口径 77 站） ----------------
# 来源: 归档 _per_rank_samples.json（同 class C、同进程数）筛选每 rank>=10 测点，
# 再用归档 cntvct_el0+none 的 met_p50_ns 剔除中位数 >=1ms 的站点；
# 运行时仍以本轮实测 met 数据复核（见 gate 阶段）。
SITES_bt="r7_l1 r7_l2 r8_l1 r9_l1 r9_l2 r9_l3 r9_l4 r10_l1 r10_l2 r10_l3 r10_l4 r11_l1 r11_l2 r11_l3 r11_l4"
SITES_cg="r3_l1 r3_l2 r3_l3 r3_l4 r3_l5 r3_l6 r3_l7 r3_l8 r4_l2"
SITES_ep="r2_l1 r3_l1"
SITES_ft="r5_l1 r6_l1 r6_l2 r7_l1 r7_l2 r7_l3 r7_l4 r15_l1"
SITES_is="r1_l2 r1_l3"
SITES_lu="r3_l1 r4_l1 r7_l1 r7_l2 r8_l1 r8_l2 r9_l1 r9_l2"
SITES_mg="r3_l1 r4_l1 r5_l1 r6_l1 r8_l1 r8_l2 r8_l3 r8_l4 r8_l5 r8_l6 r8_l7 r8_l8"
SITES_sp="r6_l1 r6_l2 r7_l1 r8_l1 r8_l2 r8_l3 r8_l4 r8_l5 r8_l6 r9_l1 r9_l2 r9_l3 r9_l4 r9_l5 r9_l6 r10_l1 r10_l2 r10_l3 r10_l4 r10_l5 r10_l6"

# 良质站点子集（--sites good）：T3（cntvct_el0+none）实测 er<0.05 且 treal 非退化
GOOD_SITES_bt="r7_l1 r7_l2 r9_l2 r9_l3 r9_l4 r10_l1 r10_l2 r10_l3 r10_l4 r11_l2 r11_l3 r11_l4"
GOOD_SITES_cg="r3_l2 r3_l5"
GOOD_SITES_ep="r2_l1"
GOOD_SITES_ft="r5_l1 r6_l1 r6_l2 r7_l1 r7_l2 r7_l3 r15_l1"
GOOD_SITES_is="r1_l2"
GOOD_SITES_lu="r3_l1 r4_l1 r7_l1 r8_l2 r9_l2"
GOOD_SITES_mg="r8_l4 r8_l5 r8_l6 r8_l7"
GOOD_SITES_sp="r6_l1 r6_l2 r8_l1 r8_l3 r8_l5 r8_l6 r9_l1 r9_l2 r9_l3 r9_l5 r9_l6 r10_l1 r10_l2 r10_l3 r10_l4 r10_l5 r10_l6"

# ---------------- 解析参数 ----------------
ONLY_BENCH=""
ONLY_SITE=""
COMBO_LIST=()
SITE_SET="full"
NO_WAIT=""
NO_FILT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY_BENCH="$2"; shift 2 ;;
    --site) ONLY_SITE="$2"; shift 2 ;;
    --combo) COMBO_LIST=($(echo "$2" | tr ',' ' ')); shift 2 ;;
    --sites) SITE_SET="$2"; shift 2 ;;
    --no-wait) NO_WAIT=1; shift ;;
    --no-filt) NO_FILT=1; shift ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done
[[ -n "$ONLY_BENCH" ]] && BENCH_LIST=($(echo "$ONLY_BENCH" | tr ',' ' '))

CONF="$SUITE/tacvar.conf"
init_root() {   # init_root TAG；TAG 非空 = 组合根（ext_<stamp>/T<tag>）
  local tag="$1"
  if [[ -n "$tag" ]]; then
    MATRIX_ROOT="$SUITE/matrix/ext_${STAMP}/${tag}"
  else
    MATRIX_ROOT="$SUITE/matrix/insitu_${STAMP}"
  fi
  DATA_ROOT_PARENT="$MATRIX_ROOT/data"
  LOG_DIR="$MATRIX_ROOT/logs"
  MANIFEST="$MATRIX_ROOT/manifest.csv"
  mkdir -p "$DATA_ROOT_PARENT" "$LOG_DIR"
  echo "$MATRIX_ROOT" > "$SUITE/matrix/latest_insitu.txt"
}

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

wait_for_idle() {  # 节点礼貌：load 超阈值则等待（--no-wait 直接放行）
  [[ -n "$NO_WAIT" ]] && return 0
  while :; do
    local ld
    ld=$(cut -d' ' -f1 /proc/loadavg)
    if awk -v l="$ld" -v m="$LOAD_MAX" 'BEGIN{exit !(l<=m)}'; then
      return 0
    fi
    log "load=${ld} > ${LOAD_MAX}, sleep ${POLL_INT}s"
    sleep "$POLL_INT"
  done
}

is_done() {  # is_done BENCH SITE ; SUCCESSFUL/MEASURED 行均视为完成
  [[ -f "$MANIFEST" ]] || return 1
  grep -qE "^${1},${2},(SUCCESSFUL|MEASURED)," "$MANIFEST"
}

check_idle_node() {  # 启动前检查：有测量进程则拒绝
  local procs
  procs=$(pgrep -f "bt.C.x|cg.C.x|ep.C.x|ft.C.x|is.C.x|lu.C.x|mg.C.x|sp.C.x|test_nspg" 2>/dev/null || true)
  if [[ -n "$procs" ]]; then
    log "REFUSE: found running NPB processes: $procs"
    exit 3
  fi
}

run_mpi() {  # run_mpi BIN [ARGS...]；NPB_TIMER_FLAG=1 + 核心绑定；stdout 由调用方重定向
  local bin="$1"; shift
  NPB_TIMER_FLAG=1 TACVAR_DATA_DIR="$data_root_abs" \
    timeout "$RUN_TIMEOUT" mpirun --map-by core --bind-to core -np "$np_this" \
    "$bin" "$@" < /dev/null
}

# ---------------- 主循环 ----------------
check_idle_node
COMBOS_RUN=("${COMBO_LIST[@]}")          # 空 => 单组合（默认行为，根=insitu_<stamp>）
[[ ${#COMBOS_RUN[@]} -eq 0 ]] && COMBOS_RUN=("")
for tag in "${COMBOS_RUN[@]}"; do
  if [[ -n "$tag" ]]; then
    IFS='|' read -r TIMER COUNTER_BACKEND COUNTER_COUNT COUNTER_NAMES NSTP <<< "${COMBO_TAGS[$tag]}"
    log "==> COMBO $tag: timer=$TIMER backend=$COUNTER_BACKEND counters=[$COUNTER_NAMES] nstp=$NSTP"
  fi
  init_root "$tag"
  if [[ ! -f "$MANIFEST" ]]; then
    echo "bench,site,status,met_verif,tfs_verif,tfe_verif,nmet,ntf,ngauge,er,data_dir,finished_at" > "$MANIFEST"
  fi
  NSPG_CACHE="$MATRIX_ROOT/nspg_cached.txt"
  NSPG_DONE=0
  log "matrix root: $MATRIX_ROOT"

for bench in "${BENCH_LIST[@]}"; do
  bench_u="${bench^^}"
  np_this="${BENCH_NP[$bench]:-$NP}"
  if [[ "$SITE_SET" == "good" ]]; then
    site_list=$(eval echo \${GOOD_SITES_${bench}})
  else
    site_list=$(eval echo \$SITES_${bench})
  fi
  [[ -z "$ONLY_SITE" ]] || site_list=$(echo "$site_list" | tr ' ' '\n' | grep -x "$ONLY_SITE" || true)
  [[ -n "$site_list" ]] || { log "no sites for $bench (filter=$ONLY_SITE), skip"; continue; }

  data_root_abs="$DATA_ROOT_PARENT/$bench"
  mkdir -p "$data_root_abs"
  bench_log="$LOG_DIR/${bench}.met.log"

  # ---- 1) 常规计时（TF=OFF）----
  if ! is_done "$bench" "$bench_u.$CLASS" || [[ ! -d "$data_root_abs/$bench.$CLASS" ]]; then
    wait_for_idle
    log "==> MET $bench.$CLASS (np=$np_this)"
    set_conf_kv TACVAR_TIMER "$TIMER"
    set_conf_kv TACVAR_NSTP "$NSTP"
    set_conf_kv TACVAR_COUNTER_BACKEND "$COUNTER_BACKEND"
    set_conf_kv TACVAR_COUNTER_COUNT "$COUNTER_COUNT"
    set_conf_kv TACVAR_COUNTER_NAMES "$COUNTER_NAMES"
    set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
    set_conf_kv TACVAR_TF_DATA_ROOT ""
    : > "$bench_log"
    met_ok=0
    if (cd "$SUITE" && make tacvar_clean && make "$bench_u" CLASS="$CLASS" PAPI_HOME="$PAPI_HOME" && make filt) >> "$bench_log" 2>&1; then
      if [[ $NSPG_DONE -eq 0 ]]; then
        if (cd "$SUITE" && make nspg NSPG_MPICC="$(command -v mpicc)" PAPI_HOME="$PAPI_HOME") >> "$bench_log" 2>&1 \
           && NPB_TIMER_FLAG=1 timeout "$RUN_TIMEOUT" mpirun --map-by core --bind-to core -np "$np_this" \
              "$SUITE/bin/test_nspg.x" "$NSPG_CACHE" >> "$bench_log" 2>&1; then
          NSPG_DONE=1
        else
          log "FAIL $bench: nspg calibration"
        fi
      fi
      if [[ $NSPG_DONE -eq 1 ]]; then
        cp "$NSPG_CACHE" "$data_root_abs/nspg.txt"
        rm -rf "$data_root_abs/$bench.$CLASS"
        if run_mpi "$SUITE/bin/${bench}.${CLASS}.x" > "$bench_log.tmp" 2>&1 \
           && grep -q SUCCESSFUL "$bench_log.tmp"; then
          cat "$bench_log.tmp" >> "$bench_log"; met_ok=1
        else
          cat "$bench_log.tmp" >> "$bench_log"
          log "FAIL $bench: met run not SUCCESSFUL"
        fi
        rm -f "$bench_log.tmp"
      fi
    else
      log "FAIL $bench: build (normal) failed"
    fi
    if [[ $met_ok -eq 1 ]]; then
      python3 "$SUITE/scripts/get_met_stat.py" "$data_root_abs" >> "$bench_log" 2>&1 \
        && echo "$bench,$bench_u.$CLASS,SUCCESSFUL,OK,,,0,0,0,,$data_root_abs,$(date '+%F %T')" >> "$MANIFEST"
    else
      echo "$bench,$bench_u.$CLASS,FAILED,,,,0,0,0,,$data_root_abs,$(date '+%F %T')" >> "$MANIFEST"
      log "skip $bench: met failed, no per-site runs"
      continue
    fi
  else
    log "skip $bench met (done)"
  fi

  # ---- 2) 站点门控：实测 per-rank 测点数 + p50 ----
  planned="$data_root_abs/planned_sites.txt"
  : > "$planned"
  for s in $site_list; do
    rid="${s#r}"; rid="${rid%_*}"
    loc="${s#*_}"; loc="${loc#l}"
    echo "$rid $loc" >> "$planned"
  done
  gate_log="$data_root_abs/gate.log"
  python3 - "$bench_u" "$CLASS" "$data_root_abs" "$planned" "$data_root_abs/sites_ok.txt" <<'PY' > "$gate_log" 2>&1
import csv, sys
import numpy as np
from collections import defaultdict
from pathlib import Path
bench_u, klass, root, planned, out = sys.argv[1:]
root = Path(root)
sites = []
with open(planned) as fp:
    for ln in fp:
        ln = ln.strip()
        if ln:
            rid, loc = ln.split()
            sites.append((int(rid), int(loc)))
vals = defaultdict(list)
per_rank = defaultdict(lambda: defaultdict(int))
met_dirs = [d for d in root.iterdir()
            if d.is_dir() and d.name.endswith("." + klass) and "_" not in d.name]
if not met_dirs:
    sys.exit(f"no met dir under {root}")
kdir = met_dirs[0]
for p in sorted(kdir.glob("*.csv")):
    if p.name == "timer_info.csv":
        continue
    with p.open(newline="") as fp:
        rd = csv.DictReader(fp)
        if not rd.fieldnames or "elapsed_ns" not in rd.fieldnames:
            continue
        for rec in rd:
            try:
                rid = int(rec["region_id"]); loc = int(rec["loc_id"])
                el = int(float(rec["elapsed_ns"]))
                rank = int(float(rec.get("rank", "0") or 0))
            except (KeyError, ValueError):
                continue
            if rec.get("migrated") in ("1", "1.0"):
                continue
            vals[(rid, loc)].append(el)
            per_rank[(rid, loc)][rank] += 1
with open(out, "w") as fp:
    for rid, loc in sites:
        arr = np.asarray(vals.get((rid, loc), []), dtype=np.int64)
        pr = per_rank.get((rid, loc), {})
        if arr.size == 0:
            print(f"r{rid}_l{loc} drop no-data", file=sys.stderr)
            continue
        minrank = min(pr.values()) if pr else 0
        p50 = int(np.median(arr))
        if minrank < 10:
            print(f"r{rid}_l{loc} drop min_per_rank={minrank}<10", file=sys.stderr)
            continue
        if p50 >= 1000000:
            print(f"r{rid}_l{loc} drop p50={p50}ns>=1ms", file=sys.stderr)
            continue
        fp.write(f"{rid} {loc} keep n={arr.size} p50={p50} minrank={minrank}\n")
        print(f"r{rid}_l{loc} keep n={arr.size} p50={p50} minrank={minrank}")
PY
  ok_sites="$data_root_abs/sites_ok.txt"
  if [[ ! -s "$ok_sites" ]]; then
    log "WARN $bench: no sites passed gate (see $gate_log)"
    continue
  fi

  # ---- 3) 逐站点：构建 TF 双二进制 + tfs/tfe 运行 + FilT ----
  # 注意: 不能用 `while read ... done < file` — 循环体里的 mpirun 会抢占 stdin 把
  # sites_ok 文件读走，导致循环只跑第一行（2026-08-15 实跑踩坑：每 kernel 只完成
  # 首站点即结束）。改为 mapfile 一次性读入列表，循环体不再共享文件 fd。
  mapfile -t sitelines < "$ok_sites"
  for ln in "${sitelines[@]}"; do
    read -r rid loc _rest <<< "$ln"
    [[ -n "$rid" ]] || continue
    site="r${rid}_l${loc}"
    if is_done "$bench" "$site"; then
      log "skip $bench $site (done)"
      continue
    fi
    wait_for_idle
    site_log="$LOG_DIR/${bench}_${site}.log"
    t0=$(date +%s)
    log "==> SITE $bench.$CLASS $site (np=$np_this)"
    set_conf_kv TACVAR_TF_SAMPLING_MODE ON
    set_conf_kv TACVAR_TF_DATA_ROOT "$data_root_abs"
    set_conf_kv TACVAR_TF_REG_ID "$rid"
    set_conf_kv TACVAR_TF_LOC_ID "$loc"
    : > "$site_log"
    ok=0
    # 3a) 构建（一次出 _tfs.x + _tfe.x）
    if (cd "$SUITE" && make tacvar_clean && make "$bench_u" CLASS="$CLASS" PAPI_HOME="$PAPI_HOME") >> "$site_log" 2>&1 \
       && [[ -x "$SUITE/bin/${bench}.${CLASS}_tfs.x" ]] && [[ -x "$SUITE/bin/${bench}.${CLASS}_tfe.x" ]]; then
      # 3b) tfs 前采样
      tfs_ok=0; tfe_ok=0
      if run_mpi "$SUITE/bin/${bench}.${CLASS}_tfs.x" > "$site_log.tmp" 2>&1 && grep -q SUCCESSFUL "$site_log.tmp"; then
        tfs_ok=1; cat "$site_log.tmp" >> "$site_log"
      else
        cat "$site_log.tmp" >> "$site_log"; log "FAIL $bench $site: tfs not SUCCESSFUL"
      fi
      # 3c) tfe 后采样
      if [[ $tfs_ok -eq 1 ]] && run_mpi "$SUITE/bin/${bench}.${CLASS}_tfe.x" > "$site_log.tmp" 2>&1 && grep -q SUCCESSFUL "$site_log.tmp"; then
        tfe_ok=1; cat "$site_log.tmp" >> "$site_log"
      else
        [[ $tfs_ok -eq 1 ]] && { cat "$site_log.tmp" >> "$site_log"; log "FAIL $bench $site: tfe not SUCCESSFUL"; }
      fi
      rm -f "$site_log.tmp"
      # 3d) FilT 过滤（run_filt.py 需 median.csv + nspg.txt 在同一 $DATA 树；
      #     --no-filt 跳过，由本机 filt_watch.sh 并行处理，站点行记 MEASURED）
      if [[ $tfs_ok -eq 1 && $tfe_ok -eq 1 ]]; then
        if [[ -n "$NO_FILT" ]]; then
          ok=1
        elif (cd "$SUITE" && python3 scripts/run_filt.py "$data_root_abs" --rid "$rid" --lid "$loc" --filt "$SUITE/bin/filt.x") >> "$site_log" 2>&1; then
          ok=1
        else
          log "FAIL $bench $site: run_filt.py"
        fi
      fi
    else
      log "FAIL $bench $site: TF build failed"
    fi

    # ---- 摘要 ----
    filt_dir="$data_root_abs/${bench}.${CLASS}_filt/r${rid}_l${loc}"
    nmet=0; ntf=0; ngauge=0; er=""
    if [[ -f "$filt_dir/met.csv" ]]; then
      nmet=$(wc -l < "$filt_dir/met.csv")
    else
      # --no-filt：直接从 tfs/tfe 目录数本站点样本（目录跨站累积，按 rid/loc 过滤）
      nmet=$(cat "$data_root_abs/${bench}.${CLASS}_tfs"/*.csv 2>/dev/null | awk -F, -v r="$rid" -v l="$loc" '$2==r && $3==l' | wc -l)
    fi
    if [[ -f "$filt_dir/tf.csv" ]]; then
      ntf=$(wc -l < "$filt_dir/tf.csv")
    else
      ntf=$(cat "$data_root_abs/${bench}.${CLASS}_tfe"/*.csv 2>/dev/null | awk -F, -v r="$rid" -v l="$loc" '$2==r && $3==l' | wc -l)
    fi
    if [[ -f "$filt_dir/er.out" ]]; then er=$(tr -d '\n' < "$filt_dir/er.out"); fi
    grep "^${bench_u},${CLASS},${rid},${loc}," "$data_root_abs/median.csv" 2>/dev/null | head -1 | awk -F, '{print $6}' > /tmp/ng.txt || true
    ngauge=$(cat /tmp/ng.txt 2>/dev/null)
    t1=$(date +%s)
    if [[ $ok -eq 1 ]]; then
      st=SUCCESSFUL; [[ -n "$NO_FILT" ]] && st=MEASURED
      echo "$bench,$site,$st,OK,,,$nmet,$ntf,${ngauge:-0},${er:-},$data_root_abs,$(date '+%F %T')" >> "$MANIFEST"
      log "OK   $bench $site [$((t1-t0))s] nmet=$nmet ntf=$ntf ngauge=${ngauge:-0} er=${er:-} (status=$st)"
    else
      echo "$bench,$site,FAILED,,,,$nmet,$ntf,${ngauge:-0},${er:-},$data_root_abs,$(date '+%F %T')" >> "$MANIFEST"
      log "FAIL $bench $site [$((t1-t0))s] see $site_log"
    fi
  done
done
done   # COMBO 循环结束

# ---- 收尾：恢复 conf 为 OFF/空，保持工作区干净 ----
set_conf_kv TACVAR_TF_SAMPLING_MODE OFF
set_conf_kv TACVAR_TF_DATA_ROOT ""
log "matrix finished. root: $MATRIX_ROOT (manifest: $MANIFEST)"
exit 0