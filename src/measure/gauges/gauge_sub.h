/**
 * @file gauge_sub.h
 * @brief Subtraction gauge and in-situ TF warmup/timed-loop macros.
 *
 * Warmup is a C register subtract (outside the timed window). The timed
 * loop continues from the warmed register for exactly nsamp decrements
 * and is emitted as architecture asm so -O2 cannot DCE it.
 */
#ifndef TACVAR_GAUGE_SUB_H
#define TACVAR_GAUGE_SUB_H

#include <stdint.h>

#if defined(__x86_64__)
#define TACVAR_TF_GAUGE_FROM_REG(ra) \
    __asm__ __volatile__( \
        "1:\n\t" \
        "sub $1, %0\n\t" \
        "jnz 1b\n\t" \
        : "+r"(ra) \
        : \
        : "cc")
#elif defined(__aarch64__)
#define TACVAR_TF_GAUGE_FROM_REG(ra) \
    __asm__ __volatile__( \
        "1:\n\t" \
        "subs %0, %0, #1\n\t" \
        "bne 1b\n\t" \
        : "+r"(ra) \
        : \
        : "cc")
#elif defined(__riscv) && (__riscv_xlen == 64)
#define TACVAR_TF_GAUGE_FROM_REG(ra) \
    __asm__ __volatile__( \
        "1:\n\t" \
        "addi %0, %0, -1\n\t" \
        "bnez %0, 1b\n\t" \
        : "+r"(ra) \
        : \
        : )
#elif defined(__loongarch64)
#define TACVAR_TF_GAUGE_FROM_REG(ra) \
    __asm__ __volatile__( \
        "1:\n\t" \
        "addi.d %0, %0, -1\n\t" \
        "bnez %0, 1b\n\t" \
        : "+r"(ra) \
        : \
        : )
#else
#define TACVAR_TF_GAUGE_FROM_REG(ra) \
    do { \
        while (ra) { \
            ra -= 1; \
        } \
    } while (0)
#endif

#define TACVAR_TF_WARMUP_REG(ra, rb, nsamp) do { \
    (ra) = (uint64_t)(nsamp) + 1ULL; \
    (rb) = 1ULL; \
    (ra) = (ra) - (rb); \
} while (0)

/** Untimed gauge of n decrements (nspg fitting). */
#define TACVAR_TF_GAUGE_N(n) do { \
    uint64_t _g_ra = (uint64_t)(n); \
    TACVAR_TF_GAUGE_FROM_REG(_g_ra); \
} while (0)

/**
 * In-situ TF sample: warmup, then time exactly nsamp decrements.
 * Requires TACVAR_TIMER_BEGIN/END from generated config.
 */
#define TACVAR_TF_SAMPLE_NS(nsamp, t0, t1) do { \
    register uint64_t _tf_ra, _tf_rb; \
    TACVAR_TF_WARMUP_REG(_tf_ra, _tf_rb, (nsamp)); \
    TACVAR_TIMER_BEGIN(t0); \
    TACVAR_TF_GAUGE_FROM_REG(_tf_ra); \
    TACVAR_TIMER_END(t1); \
    (void)_tf_rb; \
} while (0)

#endif /* TACVAR_GAUGE_SUB_H */
