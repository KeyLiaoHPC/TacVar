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
 * @param region_id NPB timer slot (or 0 for lmbench)
 * @param loc_id call-site id within the same region_id (NPB); 0 for lmbench
 */
void tacvar_csv_write_simple(int region_id, int loc_id,
                             uint64_t raw_start, uint64_t raw_stop,
                             int64_t elapsed_ns,
                             int cpu_start, int cpu_stop,
                             const uint64_t *counter_start,
                             const uint64_t *counter_stop,
                             const uint64_t *counter_delta);

/**
 * @brief Append DATA_ROOT/timer_info.csv (header once, flock around write).
 * Rows: suite,benchmark,class,test_tag,timer,region_id,nloc,name.
 * Sequential tests under the same DATA_ROOT append; concurrent writers serialize.
 * @return 0 on success, negative on error.
 */
int tacvar_write_timer_info(const int *region_ids,
                            const int *nlocs,
                            const char *const *names,
                            int n);

/** @brief Return the active data_YYYYMMDDTHHmmss directory path, or NULL. */
const char *tacvar_data_dir(void);

/** @brief True if tacvar_init completed for the current PID. */
int tacvar_is_ready(void);

/** @brief Load ngauge / TSC calibrate for in-situ TF (no-op if sampling off). */
int tacvar_tf_prepare(void);

/** @brief Subtraction-gauge loop count (round(median/nspg)) for the active site. */
uint64_t tacvar_tf_ngauge(void);

/**
 * @brief Fill per-core orig-vs-tick offset before the two gauge stamps.
 * No-op when the original timer is already TSC/CNTVCT.
 */
void tacvar_tf_ensure_offset(int cpu);

/**
 * @brief Gauge duration in ns from original-timer raw and asm-tick raw.
 * @param orig_is_start 1 for front sampling (orig then tick), 0 for rear.
 */
int64_t tacvar_tf_elapsed_ns(int orig_is_start, uint64_t t_orig, uint64_t t_tick);

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
