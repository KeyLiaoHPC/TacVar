# AGENTS.md — TacVar

Guidance for AI agents working in this repository.

## Project

TacVar measures and filters timing variability on multi-core CPUs and accelerators. Core components:

| Component | Path | Role |
|-----------|------|------|
| **ParTES** | `src/partes/` | Parallel Timing Error Sensor — Wasserstein distance between gauge runtime distributions |
| **Vkern** | `src/vkern/` | Synthetic kernels with configurable timing-fluctuation distributions |
| **Filter (FilT)** | `src/filter/` | Decouples timing fluctuation from measured runtimes |
| **Stencil** | `stencil/` | Stencil kernels and FilT integration experiments |

Related papers are listed in `README.md`. This tree is a methodology validation codebase; broader integration targets PerfHound.

## Layout

```
src/partes/     # MPI C app: gauges/, kernels/, timers/, tools/
src/filter/     # filt.c + Python helpers (get_met, get_tf, quantiles, bin width)
src/vkern/      # vkern.c, vkern_tf.c
stencil/        # stencil apps + run_*.sh
test_data/      # sample inputs
```

ParTES plugins live under `kernels/`, `gauges/`, and `timers/` and are registered via the corresponding `*.h` headers.

## Build & run

- **Language**: C99 (`-std=c99 -Wall -Wextra -O2`), plus small Python utilities under `src/filter/`.
- **ParTES**: needs MPI (`mpicc`). From `src/partes/`, set `MPI_HOME` (and `PATH` / `LD_LIBRARY_PATH` / include paths) then `make` → `partes-mpi.x`, `partes-fit.x`. Optional: `make tools`.
- **Run ParTES**: `mpirun [mpi_opts] ./partes-mpi.x --ta <ns> --tb <ns> [opts]`. See `src/partes/README.md`.
- **Filter**: build/run `filt.x` as documented in root `README.md` and `src/filter/README.md`.
- **clangd**: `.clangd` adds MPI and `src/include` paths; keep those in sync if include roots change.

Do not commit build artifacts (`.o`, `*.x`, `*.csv`, logs) — they are gitignored.

## Conventions

- Prefer existing patterns in the component you edit (ParTES Makefile dual MPI/non-MPI objects, kernel init/run/cleanup/check API, timer registration).
- New ParTES kernels/gauges/timers: implement the same function set as siblings, register in the matching header, and wire into the Makefile object lists when required.
- Keep CLI and output formats stable unless the task explicitly changes them (CSV names like `partes_*_cdf.csv`, Wasserstein/quantile tables).
- Match local C style: file-level `@file` / `@brief` comments, `pt_` / component prefixes, explicit sizes in KiB/ns where the API already uses them.
- Shell helpers (`run_*.sh`, `check_*.sh`) are the preferred way to reproduce experiments; extend them rather than inventing one-off commands when possible.

## Agent do / don't

**Do**

- Read the component README before changing CLI, build, or output.
- Build with `make` in the touched directory and smoke-test with a small `mpirun -np 2` (or the component’s script) when feasible.
- Preserve scientific meaning of metrics (W-distance, gpns, cut percentiles, bin widths).

**Don't**

- Force-push, amend shared history, or commit secrets / large result CSVs.
- Rewrite timing-critical gauge/kernel loops for “style” without measuring impact.
- Add unrelated refactors, new markdown docs, or dependency stacks unless requested.
- Assume a single global build system — each tool builds from its own directory.

## Key docs

- Root: `README.md`
- ParTES: `src/partes/README.md`
- Filter: `src/filter/README.md`
- ParTES tools: `src/partes/tools/README.md`
