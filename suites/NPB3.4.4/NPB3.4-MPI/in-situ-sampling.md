# NPB-MPI in-situ timing-fluctuation sampling

NPB-MPI is now instrumented with TacVar. When environment variable `NPB_TIMER_FLAG=1`, it can measure and save per-rank detailed timing of NPB-MPI kernels. However, these measurement results come with inevitable unstable timing fluctuation that twist the measured distribution and enlarge the heavy-tailed-end deviation. Follow this instruction to collect timing fluctuation samples and filter out timing fluctuation for more real measurement results.

Timer-adapter only. NPB kernel sources are not modified. All commands below are run from `suites/NPB3.4.4/NPB3.4-MPI` unless noted.

## 1. Flow

1. Fit nanoseconds per subtraction-gauge step (`nspg`) with `test_nspg.x`.
2. Measure the kernel with `bin/<bench>.<CLASS>.x`. Rank CSVs land under `data_<stamp>/<bench>.<CLASS>/`.
3. Copy `nspg.txt` into that data directory (the wrapper does this) and run `get_median.py`. It writes `median.csv` with `ngauge = max(10, floor(median / nspg))`.
4. Set `TACVAR_TF_SAMPLING_MODE=ON` and `TACVAR_TF_DATA_ROOT` to that directory. Rebuild: `bin/<bench>.<CLASS>_tf.x`.
5. Run the TF binary into the same data root (`TACVAR_DATA_DIR`). CSVs go to `<bench>.<CLASS>_tf/`.
6. Run `scripts/run_tf_filt.py`. For each `(region_id, loc_id)`:

   `tf = elapsed_tf - ngauge * nspg`

   FilT writes `Kernel.CLASS_filt/r{rid}_l{loc}/`, including `tf_cdf.csv` (empirical CDF of that subtracted TF).

There is one TF binary. There is no `_tfs` / `_tfe` pair and no second-pass nspg calibration.

Wrapper: `scripts/run_npb_measure.sh [measure|tf|all]` with `NP`, `BENCH`, `CLASS`.

## 2. Parameters

| Name | Default | Role |
|------|---------|------|
| `NTEST` | 10000 | Samples per `nsub` group in `test_nspg` (after warmup) |
| Warmup | 1 s `CLOCK_MONOTONIC` | Repeated `run_once(ng)` before the timed samples; not used for `tmin` |
| `ng` | 1000, 2000, 4000, … | Subtraction count per group; two `MPI_Barrier`s before each group |
| OLS window | last 5 groups | `y = a x + b` on `(nsub, tmin)`; `nspg = a` |
| `R2_THRS` | 0.999999 | Stop doubling when R² of that window reaches the threshold |
| `TMAX_NS` | 1e8 | Backstop: stop when `tmin * 2` exceeds this |
| `NSTEPS_MAX` | 32 | Backstop on number of groups |
| `ngauge` | `max(10, floor(median/nspg))` | Gauge decrements at that timer site |
| `TACVAR_TF_SAMPLING_MODE` | `OFF` | `ON` builds `*_tf.x` |
| `TACVAR_TF_DATA_ROOT` | empty | Required when `MODE=ON`; existing `data_*` directory |
| `TACVAR_TF_NSPG` | 0 | `0` means read `DATA_ROOT/nspg.txt` |
| `tf_cdf.csv` | 1000 rows | Percentiles 0.1, 0.2, …, 100.0 of subtracted TF (`percentile,tf_ns`) |

`test_nspg` starts fitting at group 5. If R² never reaches `R2_THRS`, it keeps the last valid slope and prints a warning.

## 3. How to run

Environment (fixed frequency, MPI/PAPI on `PATH`):

```bash
source scripts/setup_env.sh
```

Bind ranks to cores. Example: CG class C, 128 ranks.

**Measure and fit nspg**

```bash
# MODE=OFF in tacvar.conf
make nspg
make CG CLASS=C
mpirun --map-by core --bind-to core -np 128 ./bin/cg.C.x
mpirun --map-by core --bind-to core -np 128 ./bin/test_nspg.x data_<stamp>/nspg.txt
python3 scripts/get_median.py data_<stamp>
```

`get_median.py` requires `nspg.txt` in the data directory (or `--nspg` / `--nspg-file`).

**In-situ TF sample and FilT**

```bash
# tacvar.conf: TACVAR_TF_SAMPLING_MODE=ON
#              TACVAR_TF_DATA_ROOT=/absolute/path/to/data_<stamp>
make tacvar_clean
make CG CLASS=C
# → bin/cg.C_tf.x  and  bin/filt.x
TACVAR_DATA_DIR=data_<stamp> mpirun --map-by core --bind-to core -np 128 ./bin/cg.C_tf.x
python3 scripts/run_tf_filt.py data_<stamp> --kernel-class cg.C
```

**One-shot wrapper**

```bash
NP=128 BENCH=cg CLASS=C scripts/run_npb_measure.sh all
```

`measure` runs the kernel binary, `test_nspg`, and `get_median.py`. `tf` rebuilds `*_tf.x`, runs it, and calls FilT.

KunPeng acceptance (c920bn1, CG.C, 128 ranks, `cntvct_el0`): `bash src/measure/tests/run_arm_tests.sh --tf-c920` from the repository root.

## 4. Notes

- Use the same node, core binding, and `TACVAR_TIMER` for `test_nspg`, the kernel run, and `*_tf.x`. Mixing timers or hosts invalidates `ngauge`.
- `MODE=ON` requires `TACVAR_TF_DATA_ROOT` pointing at a directory that already contains `median.csv` (with `ngauge`) and `nspg.txt`.
- After changing `tacvar.conf`, run `make tacvar_clean` before rebuilding.
- Do not mix leftover `*_tfs.x` / `*_tfe.x` binaries or `_tfs` / `_tfe` CSV directories with this flow.
- Last `test_nspg` groups can be slow (1 s warmup plus 10000 samples at large `nsub`) if R² stays below 0.999999.
- Each FilT site directory contains `met.csv`, `tf.csv`, `tf_cdf.csv` (1000 quantile rows after subtracting `nspg * ngauge`), and FilT outputs (`tr_hist.csv`, `sim_cdf.csv`, …).
- NPB kernels under `BT/`, `CG/`, … must not be edited. Timing hooks stay in `common/c_timers.c` and `common/tacvar_npb.c`.
