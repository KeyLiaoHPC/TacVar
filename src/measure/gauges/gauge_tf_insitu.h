/**
 * @file gauge_tf_insitu.h
 * @brief Unfenced asm tick read and tick-to-ns for in-situ TF sampling.
 *
 * The inserted endpoint is always rdtscp (x86) or cntvct_el0 (Arm), with no
 * extra lfence/dmb/cpuid. Conversion matches the existing tick timers:
 * TSC * g_tacvar_ns_per_tsc, or CNTVCT * TACVAR_NSTP.
 */
#ifndef TACVAR_GAUGE_TF_INSITU_H
#define TACVAR_GAUGE_TF_INSITU_H

#include <stdint.h>
#include "gauge_sub.h"

#if defined(__x86_64__)
#include "timer_util.h"

#define TACVAR_TF_TICK_READ(raw) do { \
    unsigned _tf_hi, _tf_lo, _tf_aux; \
    __asm__ volatile("rdtscp" \
                     : "=a"(_tf_lo), "=d"(_tf_hi), "=c"(_tf_aux) \
                     :: "memory"); \
    (raw) = tacvar_join_edx_eax(_tf_hi, _tf_lo); \
} while (0)

#define TACVAR_TF_TICK_TO_NS(ticks) tacvar_tsc_delta_ns(0, (uint64_t)(ticks))

#elif defined(__aarch64__)

#ifndef TACVAR_NSTP
#define TACVAR_NSTP 0
#endif

#define TACVAR_TF_TICK_READ(raw) do { \
    uint64_t _tf_v; \
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(_tf_v) :: "memory"); \
    (raw) = _tf_v; \
} while (0)

#define TACVAR_TF_TICK_TO_NS(ticks) \
    ((int64_t)((uint64_t)(ticks) * (uint64_t)TACVAR_NSTP))

#else

#define TACVAR_TF_TICK_READ(raw) do { (raw) = 0; } while (0)
#define TACVAR_TF_TICK_TO_NS(ticks) ((int64_t)0)

#endif

#endif /* TACVAR_GAUGE_TF_INSITU_H */
