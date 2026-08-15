# DualTF 四脚本用法

把原位双采样拆成初始化 + 三阶段。四个脚本都不接收命令行参数：改 [dualtf_common.sh](dualtf_common.sh) 顶部「可编辑变量区」，然后按顺序执行。

```
bash scripts/run_dualtf_init.sh
bash scripts/run_dualtf_met.sh
bash scripts/run_dualtf_tf.sh
bash scripts/run_dualtf_filt.sh
```

可只跑其中一阶段。后三阶段在 `ROOT` 为空时读取 `matrix/latest_dualtf.txt`（由 init 写入）。

旧脚本 `run_matrix_insitu.sh` 仍保留，互不影响。

## 1 目录树

本轮根目录：`matrix/<hostname>_dualtf_YYYYMMDDTHHmmss/`（`hostname` 取变量 `NODE`）。

```
$ROOT/
  nspt.txt
  nspg.txt
  logs/
  native+none/
    met/bt.C/
    met/median.csv
    met/met_stat.csv
    tf/bt.C_r7_l1/tfs/
    tf/bt.C_r7_l1/tfe/
    filter/bt.C_r7_l1/
  native+papi_read/
  cntvct_el0+none/          # x86 自动换成 rdtscp+none
  ...
```

- `met/` 一次跑出全部站点，目录名是 `kernel.class`，没有 `_rx_lx`。
- `tf/` 与 `filter/` 按站点：`kernel.class_r${rid}_l${loc}/`。
- combo 目录名是 `<timer>+<event_reader>`，由 `TIMERS × EVENT_READERS` 生成。

## 2 变量区

共享列表在 `dualtf_common.sh` 开头。节点 MPI/PAPI/tick/默认进程数在同一文件后半，一般不用改。

| 变量 | 作用 |
|------|------|
| `NODE` | 节点名。`camd9755n2` 的 hostname 不是节点名，须改成 `camd9755n2` |
| `STAMP` / `ROOT` | init 在 `ROOT` 为空时按戳新建；后三阶段空则读 `latest_dualtf.txt` |
| `NP` | 目标进程数。`0` 用节点表（128/128/256） |
| `TIMERS` | 计时器列表，例如 `native cntvct_el0 papi_get_real_nsec` |
| `EVENT_READERS` | `none` 或 `papi_read` |
| `KERNELS` | `kernel.class` 列表，例如 `bt.C` |
| `SITES_<kernel>` | `"rid:lid"` 列表，例如 `SITES_bt=(7:1 9:2)` |
| `MPIRUN_ARGS` | 默认 `--map-by core --bind-to core -np`，后面接 `np` 和二进制 |
| `SLEEP_AFTER_KERNEL` | 每个 kernel 跑完后等待秒数，默认 5 |
| `FORCE_CALIB` | `1` 时 init 重算 nspt/nspg |
| `COMBO` | 仅 `run_dualtf_filt.sh` 顶部：要过滤的 combo 目录名 |

缩小范围（都改变量，不传参）：

- 只跑一个 combo：`TIMERS=(native)` 且 `EVENT_READERS=(none)`
- 只跑 bt.C：`KERNELS=(bt.C)`
- 只跑一站：`KERNELS=(bt.C)` 且 `SITES_bt=(9:2)`
- 换 MPI：改 `MPIRUN_ARGS` 和 `NP`

`bt` / `sp`（忽略 `.class`）用 `nearest_square_np`：取与 `NP` 距离最小的平方数，并列时取不超过 `NP` 的那个（128 → 121，256 → 256）。其它 kernel 直接用 `NP`。

x86 节点会把 `TIMERS` 里的 `cntvct_el0` 自动换成 `rdtscp`。

filt 只滤一站时，可在 `run_dualtf_filt.sh` 的 `source` 之后覆盖：

```
KERNELS=(bt.C)
SITES_bt=(7:1)
```

## 3 各阶段做什么

**init** 建目录骨架；跑 `probe_counter_freq.x --nspt-out $ROOT/nspt.txt`；用节点 tick 计时器跑一次 `test_nspg.x` 得到 `$ROOT/nspg.txt`；`make filt`。已有 nspt/nspg 且 `FORCE_CALIB=0` 则跳过标定。不跑 NPB kernel。

**met** 对每个 combo × kernel：`TF=OFF` 编译运行，不采样。`TACVAR_DATA_DIR=$combo/met`。跑完 sleep，再对 `$combo/met` 算 `median.csv` / `met_stat.csv`。已有 rank CSV 则跳过。

**tf** 对每个 combo × kernel × 站点：按 `met/median.csv` 取 ngauge，编 `_tfs.x` / `_tfe.x` 并各跑一次。写出后把 `bt.C_tfs/` 挪到 `tf/bt.C_r7_l1/tfs/`（tfe 同理）。两侧 CSV 都在则跳过。缺 median 行则该站失败并记日志。

**filt** 读 `COMBO` 目录下的 `met/` `tf/`，按 `KERNELS` × `SITES_*` 调用 `run_filt.py`，写出 `filter/kernel.class_rx_lx/`。`tf = (tfs - ngauge*nspg) + (tfe' - ngauge*nspg)`，seed=1。

## 4 注意事项

- 不要做节点 load 门控。kernel 之间只 `sleep $SLEEP_AFTER_KERNEL`。
- 站点以 `SITES_*` 为准，不再按实测 p50 或每 rank 点数丢站。
- `nspt` 和 `nspg` 只在 init 算一次。met / tf / filt 禁止重跑标定。
- `TACVAR_NSPT_FILE` 必须指向本轮 `$ROOT/nspt.txt`。缺文件时 ARM tick 会退回名义 `TACVAR_NSTP=10`。
- tfs/tfe 运行时 `TACVAR_DATA_DIR`（`$combo/tf`）里必须能读到 `median.csv` 和 `nspg.txt`，否则 `tacvar_tf_ngauge` 失败。
- 阶段 met/tf 会改 `tacvar.conf` 并 `make tacvar_clean`。不要与另一套矩阵并行编同一棵 NPB 树。
- `run_filt.py` 仍兼容旧布局（`Kernel.CLASS` 与 `_tfs/_tfe` 兄弟目录）。新树走 `met/` `tf/` `filter/`。
- 本脚本不用 git。产物 CSV 不提交。
- 在目标节点上跑。NFS 与工作台共享同一路径。
