/**
 * @file timer_util.h
 * @brief Shared helpers for TacVar timer backends (TSC calibration).
 */
#ifndef TACVAR_TIMER_UTIL_H
#define TACVAR_TIMER_UTIL_H

#ifndef _GNU_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdint.h>
#include <time.h>

#ifdef __x86_64__

/** Nanoseconds per TSC tick (set by tacvar_tsc_calibrate). */
extern double g_tacvar_ns_per_tsc;

/**
 * @brief Calibrate TSC frequency via CPUID 0x15/0x16 or CLOCK_MONOTONIC.
 * Must run before timed regions; not on the hot path.
 */
int tacvar_tsc_calibrate(void);

static inline uint64_t tacvar_join_edx_eax(unsigned hi, unsigned lo)
{
    return ((uint64_t)hi << 32) | (uint64_t)lo;
}

static inline int64_t tacvar_tsc_delta_ns(uint64_t begin, uint64_t end)
{
    uint64_t ticks = end - begin;
    return (int64_t)((double)ticks * g_tacvar_ns_per_tsc);
}

#endif /* __x86_64__ */

static inline uint64_t tacvar_timespec_ns(const struct timespec *ts)
{
    return (uint64_t)ts->tv_sec * 1000000000ULL + (uint64_t)ts->tv_nsec;
}

#endif /* TACVAR_TIMER_UTIL_H */
