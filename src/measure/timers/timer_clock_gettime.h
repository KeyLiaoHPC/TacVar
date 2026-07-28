/**
 * @file timer_clock_gettime.h
 * @brief CLOCK_MONOTONIC timer (returns nanoseconds as raw).
 */
#ifndef TACVAR_TIMER_CLOCK_GETTIME_H
#define TACVAR_TIMER_CLOCK_GETTIME_H

#ifndef _GNU_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "timer_util.h"
#include <time.h>

static inline int tacvar_timer_clock_gettime_init(void) { return 0; }
static inline void tacvar_timer_clock_gettime_fini(void) {}

static inline void tacvar_timer_clock_gettime_begin(uint64_t *raw)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    *raw = tacvar_timespec_ns(&ts);
}

static inline void tacvar_timer_clock_gettime_end(uint64_t *raw)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    *raw = tacvar_timespec_ns(&ts);
}

static inline int64_t tacvar_timer_clock_gettime_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)(e - b);
}

#endif
