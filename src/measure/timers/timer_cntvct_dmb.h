/**
 * @file timer_cntvct_dmb.h
 * @brief ARM CNTVCT_EL0 with DMB ISH barriers; raw is ticks.
 */
#ifndef TACVAR_TIMER_CNTVCT_DMB_H
#define TACVAR_TIMER_CNTVCT_DMB_H

#include <stdint.h>
#include "timer_rate.h"

#ifndef TACVAR_NSTP
#error "TACVAR_NSTP must be defined for cntvct_el0_dmb"
#endif

static inline int tacvar_timer_cntvct_dmb_init(void)
{
    return tacvar_rate_init((const char *)0, (double)TACVAR_NSTP);
}
static inline void tacvar_timer_cntvct_dmb_fini(void) {}

static inline void tacvar_timer_cntvct_dmb_begin(uint64_t *raw)
{
    uint64_t v;
    __asm__ volatile(
        "dmb ish\n\t"
        "mrs %0, cntvct_el0\n\t"
        "dmb ish\n\t"
        : "=r"(v)
        :
        : "memory");
    *raw = v;
}

static inline void tacvar_timer_cntvct_dmb_end(uint64_t *raw)
{
    uint64_t v;
    __asm__ volatile(
        "dmb ish\n\t"
        "mrs %0, cntvct_el0\n\t"
        "dmb ish\n\t"
        : "=r"(v)
        :
        : "memory");
    *raw = v;
}

static inline int64_t tacvar_timer_cntvct_dmb_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)((double)(e - b) * g_tacvar_ns_per_tick);
}

#endif
