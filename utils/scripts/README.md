# NPB-MPI timer campaign utilities

Reusable **NPB-MPI only** tooling to specify, pretest, run, summarize, and plot timer/counter comparison campaigns.

Scope: BT, CG, EP, FT, IS, LU, MG, SP. Out of scope: NPB-OMP, lmbench, agent-owned long `mpirun` jobs.

## Layout

```text
utils/scripts/
  generate_npb_timer_jobs.py      # init / pretest / fulltest generators
  run_npb_mpi_pretest.sh          # execute 00_pretest.sh
  run_npb_mpi_fulltest.sh         # generate (if needed) + execute 01_fulltest.sh
  plot_npb_mpi_step_histograms.py
  plot_npb_mpi_step_icdf.py
  analyze_npb_timer_campaign.py
  npb_timer_campaign/             # importable package
  pyproject.toml / uv.lock / environment.yaml
```

Skill (agent workflow): [`SKILLS/test-npb-mpi-timer-comparison/SKILL.md`](../../SKILLS/test-npb-mpi-timer-comparison/SKILL.md).

## Environment setup

```bash
cd utils/scripts
uv sync
source .venv/bin/activate
```

[`environment.yaml`](environment.yaml) lists the same Python ≥3.10 + matplotlib constraint for humans who do not use uv. After sync, plain `python3` works inside `.venv`.

Publication figures require **Times New Roman**. If matplotlib cannot resolve it, plotting fails unless you pass `--allow-font-fallback` / `allow_font_fallback=True` (recorded in `REPORT.md`).

## Required campaign fields

No silent defaults. Collect explicitly:

| Field | Notes |
|---|---|
| `hostname` | Target login/compute host |
| `cores` | Default `mpirun -np` |
| `kernels` | Subset of BT,CG,EP,FT,IS,LU,MG,SP |
| `timers` | Valid for host arch |
| `counter_profiles` | At least one `{name,backend,events[]}`; e0 uses `backend=none`, `events=[]` |
| `class` | NPB class |
| `mpi_home`, `papi_home` | Absolute paths on target host |
| `bind_args` | e.g. `--map-by core --bind-to core` |
| `kernel_np` | Required for BT/SP when `cores` is not square |

## Two-stage workflow

### 1. Freeze specification + generate pretest

```bash
cd utils/scripts
python3 generate_npb_timer_jobs.py init \
  --hostname c920bn1 \
  --cores 128 \
  --kernels BT,CG,EP,FT,IS,LU,MG,SP \
  --timers native,clock_gettime,papi_get_real_nsec,cntvct_el0,cntvct_el0_dmb \
  --profiles-json '[{"name":"e0","backend":"none","events":[]},{"name":"e4","backend":"papi_read","events":["PAPI_TOT_CYC","PAPI_TOT_INS","PAPI_L1_DCM","PAPI_LST_INS"]}]' \
  --class C \
  --mpi-home /home/hpckey/01-App/openmpi-5.0.8 \
  --papi-home /home/hpckey/01-App/papi \
  --square-np 121 \
  --arch aarch64
```

Or from an existing `campaign.json`:

```bash
python3 generate_npb_timer_jobs.py pretest --spec /path/to/campaign.json
```

### 2. Run pretest (detached; user-owned)

```bash
ssh <host> "cd <results-dir> && nohup bash 00_pretest.sh > pretest.nohup.out 2>&1 &"
# or:
./run_npb_mpi_pretest.sh <results-dir>
```

Expect `preflight/STATUS=PASS` and `preflight/spec.sha256` matching `campaign.sha256`. Failed events stay failed (no hidden fallback).

### 3. Generate + run fulltest only after matching PASS

```bash
python3 generate_npb_timer_jobs.py fulltest --results-dir <results-dir>
ssh <host> "cd <results-dir> && nohup bash 01_fulltest.sh > fulltest.nohup.out 2>&1 &"
# or:
./run_npb_mpi_fulltest.sh --results-dir <results-dir> --generate
```

Resumable: cells with `STATUS=PASS` and matching `run-manifest.json` digest are skipped. Scripts restore `tacvar.conf` via EXIT trap and use `/tmp/tacvar-npb-mpi-timer-${USER}.lock`.

## Analysis and plotting

Supports:

1. New trees with `campaign.json`
2. Legacy trees such as  
   `suites/NPB3.4.4/NPB3.4-MPI/results_c920bn1_classC_20260728T181952/`  
   (matrix inferred from `runs/<timer>_<profile>/<kernel>/`)

Per-step samples are TacVar CSV rows with `region_id=1000`.

```bash
RESULTS=suites/NPB3.4.4/NPB3.4-MPI/results_c920bn1_classC_20260728T181952

python3 plot_npb_mpi_step_histograms.py "$RESULTS"
python3 plot_npb_mpi_step_icdf.py "$RESULTS"
python3 analyze_npb_timer_campaign.py "$RESULTS"   # summaries + both plots + REPORT.md
```

Filters: `--kernels`, `--timers`, `--profiles`. Incomplete matrix refuses analysis unless `--allow-partial`.

Outputs under `summary/`:

- `overall_runtime.csv`, `total_dist_stats.csv`, `detail_from_stdout.csv`, `step_dist_stats.csv`
- `step_histograms.json`, `step_icdf.json`, `step_dist_analysis.json`
- `figures/hist_<kernel>_<profile>.{png,pdf,eps}`
- `figures/icdf_<kernel>_<profile>.{png,pdf,eps}`
- `REPORT.md` (relative PNG embeds; PDF/EPS links)

### Notebook / import API

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path("/path/to/TacVar/utils/scripts")))

from npb_timer_campaign import (
    analyze_campaign,
    load_campaign,
    plot_step_histograms,
    plot_step_icdf,
    summarize_campaign,
    write_report,
)

results = Path("/path/to/results_...")
campaign = load_campaign(results)
summary = summarize_campaign(campaign)
hist = plot_step_histograms(campaign, formats=("png", "pdf", "eps"))
icdf = plot_step_icdf(campaign, formats=("png", "pdf", "eps"))
write_report(campaign, summary, hist, icdf)
# or:
analyze_campaign(results, formats=("png", "pdf", "eps"))
```

## Example: current Class C fixture (not a default)

Request that produced the reference tree:

- Host: `c920bn1` (aarch64)
- Class C, `mpirun -np 128 --map-by core --bind-to core` (BT/SP → 121)
- Timers: `native`, `clock_gettime`, `papi_get_real_nsec`, `cntvct_el0`, `cntvct_el0_dmb`
- Profiles: `e0` (`none`) and `e4` (`papi_read` with `PAPI_TOT_CYC,PAPI_TOT_INS,PAPI_L1_DCM,PAPI_LST_INS`)
- Results: `suites/NPB3.4.4/NPB3.4-MPI/results_c920bn1_classC_20260728T181952/`

Compatibility wrappers (labeled legacy defaults):

- `suites/NPB3.4.4/run_npb_timer_classC_c920.sh` → analyze or legacy one-shot
- `suites/NPB3.4.4/summarize_npb_timer_classC.py` → `analyze_npb_timer_campaign.py`

## Troubleshooting

| Symptom | Action |
|---|---|
| Times New Roman missing | Install msttcorefonts / supply TNR to matplotlib, or `--allow-font-fallback` |
| fulltest generator refuses | Ensure `preflight/STATUS=PASS` and digests match |
| Zero event deltas | Fix events in `campaign.json` and rerun pretest (no silent substitution) |
| Incomplete matrix | Finish FAIL cells or use `--allow-partial` for diagnostics |
| Lock busy / conflicting mpirun | Inspect reported PIDs; scripts never auto-kill |

## Data retention / gitignore

Generated `results_*`, `.venv`, Python/plot caches, and `tacvar.conf.*.bak` are gitignored. Tracked: `pyproject.toml`, `uv.lock`, `environment.yaml`, package sources, this README, and the skill.
