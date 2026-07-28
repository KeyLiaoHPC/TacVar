/**
 * @file timer_rdtsc.h
 * @brief Non-serializing RDTSC; raw is TSC ticks.
 */
#ifndef TACVAR_TIMER_RDTSC_H
#define TACVAR_TIMER_RDTSC_H

#include "timer_util.h"

static inline int tacvar_timer_rdtsc_init(void)
{
    return tacvar_tsc_calibrate();
}
static inline void tacvar_timer_rdtsc_fini(void) {}

static inline void tacvar_timer_rdtsc_begin(uint64_t *raw)
{
    unsigned hi, lo;
    __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi) :: "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline void tacvar_timer_rdtsc_end(uint64_t *raw)
{
    unsigned hi, lo;
    __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi) :: "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline int64_t tacvar_timer_rdtsc_delta_ns(uint64_t b, uint64_t e)
{
    return tacvar_tsc_delta_ns(b, e);
}

#endif
