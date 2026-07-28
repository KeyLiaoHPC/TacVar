/**
 * @file timer_rdtscp_lfence.h
 * @brief LFENCE+RDTSC begin / RDTSCP+LFENCE end; raw is TSC ticks.
 */
#ifndef TACVAR_TIMER_RDTSCP_LFENCE_H
#define TACVAR_TIMER_RDTSCP_LFENCE_H

#include "timer_util.h"

static inline int tacvar_timer_rdtscp_lfence_init(void)
{
    return tacvar_tsc_calibrate();
}
static inline void tacvar_timer_rdtscp_lfence_fini(void) {}

static inline void tacvar_timer_rdtscp_lfence_begin(uint64_t *raw)
{
    unsigned hi, lo;
    __asm__ volatile("lfence\n\trdtsc"
                     : "=a"(lo), "=d"(hi)
                     :
                     : "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline void tacvar_timer_rdtscp_lfence_end(uint64_t *raw)
{
    unsigned hi, lo, aux;
    __asm__ volatile("rdtscp\n\tlfence"
                     : "=a"(lo), "=d"(hi), "=c"(aux)
                     :
                     : "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline int64_t tacvar_timer_rdtscp_lfence_delta_ns(uint64_t b, uint64_t e)
{
    return tacvar_tsc_delta_ns(b, e);
}

#endif
