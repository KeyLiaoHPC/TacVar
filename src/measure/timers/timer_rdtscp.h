/**
 * @file timer_rdtscp.h
 * @brief RDTSCP (serializes prior instructions); raw is TSC ticks.
 */
#ifndef TACVAR_TIMER_RDTSCP_H
#define TACVAR_TIMER_RDTSCP_H

#include "timer_util.h"

static inline int tacvar_timer_rdtscp_init(void)
{
    return tacvar_tsc_calibrate();
}
static inline void tacvar_timer_rdtscp_fini(void) {}

static inline void tacvar_timer_rdtscp_begin(uint64_t *raw)
{
    unsigned hi, lo, aux;
    __asm__ volatile("rdtscp" : "=a"(lo), "=d"(hi), "=c"(aux) :: "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline void tacvar_timer_rdtscp_end(uint64_t *raw)
{
    unsigned hi, lo, aux;
    __asm__ volatile("rdtscp" : "=a"(lo), "=d"(hi), "=c"(aux) :: "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline int64_t tacvar_timer_rdtscp_delta_ns(uint64_t b, uint64_t e)
{
    return tacvar_tsc_delta_ns(b, e);
}

#endif
