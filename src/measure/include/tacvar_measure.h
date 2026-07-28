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
 *
 * Event rows are buffered in memory and written only during tacvar_fini().
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

/** Reserved region_id for NPB major-loop / per-step intervals. */
#ifndef TACVAR_REGION_STEP
#define TACVAR_REGION_STEP 1000
#endif

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

/** Logical region metadata for region_info.csv. */
typedef struct tacvar_region_info {
    int region_id;
    const char *region_name;
    const char *source_locations; /**< "file:routine:start-end;..." */
    const char *description;
    const char *active_when;
    int recorded_by_tacvar;       /**< 1 if written to event CSV */
} tacvar_region_info_t;

/**
 * @brief Initialize measurement (directory, backends; path only for CSV).
 * Safe to call multiple times in the same process; re-inits after fork.
 * @return 0 on success, negative errno-style code on failure.
 */
int tacvar_init(const tacvar_context_t *ctx);

/** @brief Flush buffered events / metadata and release backends. */
void tacvar_fini(void);

/**
 * @brief Append one measured interval to the in-memory event buffer.
 * Must be called after TACVAR_TIMER_END (outside the timed region).
 * Rows are written to disk only in tacvar_fini().
 * @param source  Caller source location string, or NULL.
 */
void tacvar_csv_write_simple(int region_id,
                             uint64_t raw_start, uint64_t raw_stop,
                             int64_t elapsed_ns,
                             int cpu_start, int cpu_stop,
                             const uint64_t *counter_start,
                             const uint64_t *counter_stop,
                             const uint64_t *counter_delta,
                             const char *source);

/**
 * @brief Preallocate capacity for exactly required_events.
 * @return 0 on success, negative errno-style code on failure.
 */
int tacvar_events_reserve(size_t required_events);

/** @brief Clear registered region metadata for the current writer. */
void tacvar_region_info_clear(void);

/** @brief Register logical regions (copied by value into internal table). */
void tacvar_region_info_register(const tacvar_region_info_t *regions,
                                 size_t count);

/**
 * @brief Write <data_dir>/region_info.csv once (rank/thread 0).
 * Idempotent if the file already exists.
 * @return 0 on success / skip, negative on failure.
 */
int tacvar_region_info_write(void);

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
