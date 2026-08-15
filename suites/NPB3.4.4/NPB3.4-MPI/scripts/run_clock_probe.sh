#!/usr/bin/env bash
# run_clock_probe.sh — c920bn1 上执行 79ppm 时钟偏差稳定性测量
#   1) gauge 级短间隔: 多个 ng 的减法 gauge 双时钟比率 + 分轮跨时间稳定序列
#   2) 长间隔: 真实 128 进程 NPB(lu.C) 运行前后读两时钟绝对值, 计算全程比率与 offset 漂移
set -uo pipefail
SUITE=/mnt/keylabmain/nfs/hpckey/03-Project/TacVar_NPC/suites/NPB3.4.4/NPB3.4-MPI
cd "$SUITE" || exit 2
export MPI_HOME=/home/hpckey/01-App/openmpi-5.0.8
export PAPI_HOME=/home/hpckey/01-App/papi
export PATH="$MPI_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_HOME/lib:$PAPI_HOME/lib:${LD_LIBRARY_PATH:-}"
PROBE=bin/probe_clock_rate
OUT=/tmp/clockprobe
mkdir -p "$OUT"

echo "=== build probe ==="
cc -O2 -o "$PROBE" scripts/probe_clock_rate.c -lm || exit 3

echo
echo "=== [1] gauge-based: 多 ng 短间隔比率 ==="
"$PROBE" --gauge 51304 2000
"$PROBE" --gauge 513040 500
"$PROBE" --gauge 3000000 60

echo
echo "=== [2] gauge-based: 跨时间稳定序列 (ng=513040, 每轮 300 个, 轮间 10s) ==="
"$PROBE" --series 513040 300 6 10

echo
echo "=== [3] build LU CLASS=C (TF=OFF) ==="
make tacvar_clean >/dev/null 2>&1
make LU CLASS=C PAPI_HOME="$PAPI_HOME" > "$OUT/lu_build.log" 2>&1 || { tail -20 "$OUT/lu_build.log"; exit 4; }
echo "build ok"

echo
echo "=== [4] long-interval: 128-rank lu.C 前后绝对读数 ==="
"$PROBE" --abs > "$OUT/before.txt"
cat "$OUT/before.txt"
echo "--- run lu.C (np=128) ---"
NPB_TIMER_FLAG=1 mpirun --map-by core --bind-to core -np 128 ./bin/lu.C.x > "$OUT/lu_run.log" 2>&1
grep -E "Verification|Time in seconds" "$OUT/lu_run.log" | head -3
"$PROBE" --abs > "$OUT/after.txt"
echo "--- after ---"
cat "$OUT/after.txt"

echo
echo "=== [5] 计算 ==="
python3 - "$OUT/before.txt" "$OUT/after.txt" <<'PY'
import sys
import numpy as np
def load(p):
    rows = []
    with open(p) as fp:
        for ln in fp:
            g, c = ln.split()
            rows.append((int(g), int(c)))
    return np.array(rows, dtype=np.int64)
b = load(sys.argv[1]); a = load(sys.argv[2])
# 用中位绝对对（排序后配对，减少读序噪声）
bg, bc = np.median(b[:,0]), np.median(b[:,1])
ag, ac = np.median(a[:,0]), np.median(a[:,1])
d_g = ag - bg; d_c = ac - bc
ratio = (d_c * 10.0) / d_g
gap_b = bg - bc * 10
gap_a = ag - ac * 10
print(f"long-interval: elapsed={d_g/1e9:.3f} s (clock); {d_c} cntvct ticks")
print(f"  ratio (cntvct-ns/clock-ns) = {ratio:.6f}  => 偏差 {(1-ratio)*1e6:.1f} ppm")
print(f"  offset gap before={gap_b/1e9:.6f} s  after={gap_a/1e9:.6f} s")
print(f"  gap drift over run = {(gap_a-gap_b)/1e9:.6f} s = {(gap_a-gap_b)/d_g*1e6:.1f} ppm of elapsed")
PY