---
name: test-npb-mpi-timer-comparison
description: >-
  Two-stage NPB-MPI timer/counter comparison campaigns: collect an explicit
  specification, generate and hand off a guarded pretest, generate fulltest only
  after matching PASS, then summarize, plot, and write REPORT.md. Use when the
  user asks to compare NPB-MPI timers, counters, Class C/A/B runs, or per-step
  timing distributions across timers.
---

# NPB-MPI timer comparison skill

## When to use

User wants an NPB-MPI experiment comparing timers and/or counter profiles (total, NPB detail, TacVar per-step), including analysis of existing `results_*` trees.

## Mandatory reads

1. [`utils/scripts/README.md`](../../utils/scripts/README.md)
2. This skill

## Hard rules

1. **Ask for missing fields.** Never silently choose hostname, cores, kernels, timers, counter readers/events, class, MPI/PAPI paths, or bind args.
2. **NPB-MPI only** in this skill (no OMP/lmbench).
3. **Two-stage generation.** Do not generate `01_fulltest.sh` until `preflight/STATUS=PASS` and `preflight/spec.sha256` matches `campaign.sha256`.
4. **No agent-owned long mpirun.** Generate scripts and return detached launch commands. Do not keep matrix jobs attached to the agent tool lifecycle unless the user explicitly asks to monitor.
5. **No hidden event fallbacks.** Failed requested events remain failed; changing events means updating the spec and rerunning pretest.
6. **Never auto-kill** conflicting processes; report them.

## Workflow

### Stage A — collect + freeze

Ask only for missing values from:

- hostname, cores, kernels, timers, counter_profiles, class
- MPI_HOME, PAPI_HOME, bind_args
- BT/SP square `kernel_np` when cores are not square

Write `campaign.json` / digest and generate pretest:

```bash
cd utils/scripts && uv sync
python3 generate_npb_timer_jobs.py init ...   # or: pretest --spec campaign.json
```

Return the pretest launch command and **stop**:

```bash
ssh <host> "cd <results-dir> && nohup bash 00_pretest.sh > pretest.nohup.out 2>&1 &"
```

### Stage B — after user confirms pretest finished

Verify PASS + digest, then:

```bash
python3 utils/scripts/generate_npb_timer_jobs.py fulltest --results-dir <results-dir>
```

Return the fulltest launch command and **stop**.

### Stage C — after FULLTEST_STATUS=COMPLETE

```bash
python3 utils/scripts/analyze_npb_timer_campaign.py <results-dir>
```

Always produce `REPORT.md` + PNG/PDF/EPS figures. Create a Cursor canvas only when canvas support is available; otherwise say so and keep markdown/figures as the archive.

## Analysis of existing trees

Legacy fixture (no `campaign.json`):

`suites/NPB3.4.4/NPB3.4-MPI/results_c920bn1_classC_20260728T181952/`

Discovery infers kernels/timers/profiles from `runs/`. Do not hardcode those lists in new plot code.

## Prohibitions

- Ad-hoc one-off matrix shell that bypasses `campaign.json`
- Inventing Class C / 128 / c920 defaults when the user did not state them
- Committing `.venv`, CSVs, or `results_*`
