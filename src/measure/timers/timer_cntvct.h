/**
 * @file timer_cntvct.h
 * @brief ARM CNTVCT_EL0; raw is ticks; ns = ticks * TACVAR_NSTP.
 */
#ifndef TACVAR_TIMER_CNTVCT_H
#define TACVAR_TIMER_CNTVCT_H

#include <stdint.h>

#ifndef TACVAR_NSTP
#error "TACVAR_NSTP must be defined for cntvct_el0"
#endif

static inline int tacvar_timer_cntvct_init(void) { return 0; }
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
    return (int64_t)((e - b) * (uint64_t)TACVAR_NSTP);
}

#endif
