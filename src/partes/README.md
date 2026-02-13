# ParTES - Parallel Timing Error Sensor

## 1 Overview

ParTES (Parallel Timing Error Sensor) measures the deviation between two measured runtime distributions. It targets a common but ambiguous practice in parallel performance measurement: enlarging the timing interval to reduce relative timing error. This can fail for two reasons: 1) adding instructions to enlarge the interval distorts per-step runtime distributions and can introduce new variability; and 2) there is no principled way to decide how large is “large enough” to make errors negligible. ParTES provides a systematic way to quantify timing error across runtime scenarios and timing intervals, helping users make informed measurement decisions.

In ParTES, each distribution represents aggregated per-core, per-step measurements of a “gauge.” A gauge is a loop kernel that repeatedly executes one instruction (e.g., subtraction) and is designed to avoid disturbing memory or the network subsystem. The goal is to make the gauge the most stable kernel with minimal variability. By adding other operations before and/or after the gauge, ParTES quantifies timing deviation in “noisy” run-time environments with cache misses, frequency changes, process preemption, and similar effects.

ParTES uses the 1-D Wasserstein distance to measure differences between two runtime distributions, ta and tb. Users set run-time parameters through the ParTES CLI (ta, tb, kernel and flush sizes, etc.) to quantify per-core timing error across timing intervals, run-time scenarios, timers, parallel scales, and surrounding kernels.

## 2 Structure

### 2.1 Kernel Files

All kernels are located in the `kernels/` directory:

- `none.c` - No operation kernel
- `triad.c` - Triad operation: `a[i] = 0.42 * b[i] + c[i]`
- `scale.c` - Scale operation: `a[i] = 1.0001 * b[i]; b[i] = 1.0001 * a[i]`
- `copy.c` - Copy operation: `a[i] = b[i]`
- `add.c` - Add operation: `a[i] = b[i] + c[i]`
- `pow.c` - Power operation: `a[i] = pow(b[i], 1.0001)`
- `dgemm.c` - Matrix multiplication (simplified)
- `bcast.c` - MPI broadcast operation

### 2.2 Each Kernel Implements

1. **Init function** - Initialize kernel data structures with memory size in KiB
2. **Run function** - Execute the kernel operation
3. **Cleanup function** - Free allocated memory

## 3 Usage

### 3.1 Build

**Prerequisites**
- MPI library (e.g., MPICH, OpenMPI, MVAPICH, Intel MPI)
- C compiler with MPI support (e.g., gcc >= 8.5)

**Environmental Variables**

Assume the MPI library is installed in `/opt/mpi`:
```bash
export MPI_HOME=/opt/mpi
export PATH=$MPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$MPI_HOME/lib:$LIBRARY_PATH
export CPATH=$MPI_HOME/include:$CPATH
export C_INCLUDE_PATH=$MPI_HOME/include:$C_INCLUDE_PATH
```
Then type `make` to build the project. After compilation, the executable `partes-mpi.x` will be generated in the current directory.

### 3.2 Run `partes-mpi.x`

Use `mpirun`, `mpiexec` or related MPI launcher to run `partes-mpi.x`.
```bash
mpirun [mpi_options] ./partes-mpi.x [partes_options]
```

Mandatory partes options:
- `--ta <time(ns)>`: Theoretical execution time in ns of the first gauge.
- `--tb <time(ns)>`: Theoretical execution time in ns of the second gauge.

Optional partes options:
- `--fkern-a <kernel>`: Front kernel, default: none. Refer to `./kernels/kernels.h` for available kernels (none, triad, scale, copy, add, pow, dgemm, mpi_bcast).
- `--fsize-a <size>`: Front kernel memory size in KiB, default: 0.
- `--fkern-b <kernel>`: Rear kernel, default: none.
- `--fsize-b <size>`: Rear kernel memory size in KiB, default: 0.
- `--rkern-a <kernel>`: Front kernel for second gauge, default: none.
- `--rsize-a <size>`: Rear kernel memory size for second gauge in KiB, default: 0.
- `--rkern-b <kernel>`: Rear kernel for second gauge, default: none.
- `--rsize-b <size>`: Rear kernel memory size for second gauge in KiB, default: 0.
- `--timer <timer>`: Timer method (clock_gettime, mpi_wtime)
- `--ntests <num>`: Number of gauge measurements (default: 1000)
- `--timer <timing_method>`: Timing method (default: clock_gettime), refer to `timers/timers.h` for available timers
- `--gauge <gauge_kernel>`: Gauge kernel (default: sub_scalar), refer to `gauges/gauges.h` for available gauges.
- `--ntiles <num>`: Number of tiles (default: 100).
- `--cut-p <num>`: Percentage cut for outlier removal (default: 1.0).
- `--help, -h`: Show help message

### 3.3 Outputs and examples

**Example 1**:
Run tests to measure the W-Distance of clock_gettime timing functions under the 2-core per-core timing measurement scenario. Set theoretical of two gauges to 1000ns and 2000ns, without flush kernels. Here in the outputs, the line leading by "#---" is the comment line, which is not printed in the actual output.
```bash
$ mpirun -np 2 ./partes-mpi.x --ta 1000 --tb 2000
#--- Run-time settings
Repeat 1000 runtime measurements, target gauge time: 1000ns, 2000ns
Timer: clock_gettime
Gauge: sub_scalar
ta flush info:
Front kernel: NONE, size: 0 KiB, real size: 0 KiB
Rear kernel: NONE, size: 0 KiB, real size: 0 KiB
tb flush info:
Front kernel: NONE, size: 0 KiB, real size: 0 KiB
Rear kernel: NONE, size: 0 KiB, real size: 0 KiB
#--- Considering the main frequency of the processor could be volatile, ParTES measures the number of gauge step per nanosecond (gpns).
[Rank 0] Gauge info: gpns=3.098000
[Rank 1] Gauge info: gpns=3.098000
t0 = 1000, number of gauges: 3098
t1 = 2000, number of gauges: 6196
#--- Result checks for avoiding the overly optimization of the flush kernel.
TA Front kernel percentage gap: 0.000000%
TA Rear kernel percentage gap: 0.000000%
TB Front kernel percentage gap: 0.000000%
TB Rear kernel percentage gap: 0.000000%
#--- The tailed fraction being cutted.
Percentage cut: 1.000000
#--- Theoretical time gap between ta and tb
Time gap: 1000ns
#--- Measured W per key quantile.
Quantile, W(Ta), W(Tb), W(Tb)-W(Ta)
0, 1027, 2026, 999
50, 1030, 2029, 999
75, 1031, 2030, 999
90, 1031, 2031, 1000
95, 1032, 2040, 1008
99, 2109, 887338, 885229
100, 2109, 887338, 885229
#--- Overall W
Wasserstein distance: 9842.330000
```
The per-step measured running time of each core is saved to `pates_<ta/tb>_r<rank_id>.csv`. The cumulative density function discrete array of the first and second gauges' measured running times is saved to `partes_<ta/tb>_cdf.csv` in the current directory.

