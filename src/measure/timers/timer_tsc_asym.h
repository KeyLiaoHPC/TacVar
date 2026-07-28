/**
 * @file timer_tsc_asym.h
 * @brief Paoloni asymmetric TSC: CPUID+RDTSC / RDTSCP+CPUID.
 * Raw is TSC ticks. Correct packing: ((uint64_t)edx << 32) | eax.
 */
#ifndef TACVAR_TIMER_TSC_ASYM_H
#define TACVAR_TIMER_TSC_ASYM_H

#include "timer_util.h"

static inline int tacvar_timer_tsc_asym_init(void)
{
    return tacvar_tsc_calibrate();
}
static inline void tacvar_timer_tsc_asym_fini(void) {}

static inline void tacvar_timer_tsc_asym_begin(uint64_t *raw)
{
    unsigned hi, lo;
    __asm__ volatile(
        "cpuid\n\t"
        "rdtsc\n\t"
        : "=a"(lo), "=d"(hi)
        : "a"(0)
        : "rbx", "rcx", "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
}

static inline void tacvar_timer_tsc_asym_end(uint64_t *raw)
{
    unsigned hi, lo, aux;
    __asm__ volatile(
        "rdtscp\n\t"
        : "=a"(lo), "=d"(hi), "=c"(aux)
        :
        : "memory");
    *raw = tacvar_join_edx_eax(hi, lo);
    __asm__ volatile(
        "cpuid\n\t"
        :
        : "a"(0)
        : "rbx", "rcx", "rdx", "memory");
}

static inline int64_t tacvar_timer_tsc_asym_delta_ns(uint64_t b, uint64_t e)
{
    return tacvar_tsc_delta_ns(b, e);
}

#endif
