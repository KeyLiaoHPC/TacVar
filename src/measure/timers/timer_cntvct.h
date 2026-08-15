/**
 * @file timer_cntvct.h
 * @brief ARM CNTVCT_EL0; raw is ticks; ns = ticks * g_tacvar_ns_per_tick.
 */
#ifndef TACVAR_TIMER_CNTVCT_H
#define TACVAR_TIMER_CNTVCT_H

#include <stdint.h>
#include "timer_rate.h"

#ifndef TACVAR_NSTP
#error "TACVAR_NSTP must be defined for cntvct_el0"
#endif

static inline int tacvar_timer_cntvct_init(void)
{
    return tacvar_rate_init((const char *)0, (double)TACVAR_NSTP);
}
static inline void tacvar_timer_cntvct_fini(void) {}

static inline void tacvar_timer_cntvct_begin(uint64_t *raw)
{
    uint64_t v;
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(v) :: "memory");
    *raw = v;
}

static inline void tacvar_timer_cntvct_end(uint64_t *raw)
{
    uint64_t v;
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(v) :: "memory");
    *raw = v;
}

static inline int64_t tacvar_timer_cntvct_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)((double)(e - b) * g_tacvar_ns_per_tick);
}

#endif
