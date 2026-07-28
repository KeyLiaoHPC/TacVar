/**
 * @file tacvar_measure.h
 * @brief Public TacVar measurement API for NPB and lmbench.
 *
 * Hot-path timer/counter reads use compile-time macros generated from
 * tacvar.conf (TACVAR_TIMER_BEGIN/END/DELTA_NS, TACVAR_COUNTER_READ).
 * This header only exposes non-hot-path init/fini and CSV helpers.
 *
 * Units: timer raw values are backend-specific; TACVAR_TIMER_DELTA_NS
 * always returns nanoseconds. Counter values are raw hardware counts.
 */
#ifndef TACVAR_MEASURE_H
#define TACVAR_MEASURE_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum supported hardware counters in one build. */
#define TACVAR_MAX_COUNTERS 12

/** Context passed to tacvar_init(); all strings may be NULL. */
typedef struct tacvar_context {
    const char *suite;       /**< "npb-mpi", "npb-omp", "lmbench" */
    const char *benchmark;   /**< e.g. "is", "cg", "lat_syscall" */
    const char *klass;       /**< NPB class letter, or NULL */
    const char *test_tag;    /**< optional sub-test tag */
    int rank;                /**< MPI rank, else 0 */
    int thread;              /**< OpenMP thread id, else 0 */
    int nprocs;              /**< MPI size / OMP threads, else 1 */
} tacvar_context_t;

/**
 * @brief Initialize measurement (directory, CSV, selected backends).
 * Safe to call multiple times in the same process; re-inits after fork.
 * @return 0 on success, negative errno-style code on failure.
 */
int tacvar_init(const tacvar_context_t *ctx);

/** @brief Flush CSV and release backend resources for the current PID. */
void tacvar_fini(void);

/**
 * @brief Append one measured interval to the per-writer CSV.
 * Must be called after TACVAR_TIMER_END (outside the timed region).
 */
void tacvar_csv_write_simple(int region_id,
                             uint64_t raw_start, uint64_t raw_stop,
                             int64_t elapsed_ns,
                             int cpu_start, int cpu_stop,
                             const uint64_t *counter_start,
                             const uint64_t *counter_stop,
                             const uint64_t *counter_delta);

/** @brief Return the active data_YYYYMMDDTHHmmss directory path, or NULL. */
const char *tacvar_data_dir(void);

/** @brief True if tacvar_init completed for the current PID. */
int tacvar_is_ready(void);

#ifdef __cplusplus
}
#endif

/* Generated config supplies hot-path macros when present. */
#if defined(TACVAR_GENERATED_CONFIG_H)
#include TACVAR_GENERATED_CONFIG_H
#elif defined(__has_include)
#  if __has_include("tacvar_generated_config.h")
#    include "tacvar_generated_config.h"
#  endif
#else
/* Fallback: try relative include used by suite build dirs. */
#endif

#endif /* TACVAR_MEASURE_H */
