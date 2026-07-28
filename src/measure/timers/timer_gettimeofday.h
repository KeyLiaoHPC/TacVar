/**
 * @file timer_gettimeofday.h
 * @brief gettimeofday timer (lmbench native); raw is nanoseconds.
 */
#ifndef TACVAR_TIMER_GETTIMEOFDAY_H
#define TACVAR_TIMER_GETTIMEOFDAY_H

#include <sys/time.h>
#include <stdint.h>

static inline int tacvar_timer_gettimeofday_init(void) { return 0; }
static inline void tacvar_timer_gettimeofday_fini(void) {}

static inline void tacvar_timer_gettimeofday_begin(uint64_t *raw)
{
    struct timeval tv;
    gettimeofday(&tv, (struct timezone *)0);
    *raw = (uint64_t)tv.tv_sec * 1000000000ULL
         + (uint64_t)tv.tv_usec * 1000ULL;
}

static inline void tacvar_timer_gettimeofday_end(uint64_t *raw)
{
    struct timeval tv;
    gettimeofday(&tv, (struct timezone *)0);
    *raw = (uint64_t)tv.tv_sec * 1000000000ULL
         + (uint64_t)tv.tv_usec * 1000ULL;
}

static inline int64_t tacvar_timer_gettimeofday_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)(e - b);
}

#endif
