/**
 * @file timer_rate.h
 * @brief Platform-independent nanosecond-per-tick rate for cycle counters.
 *
 * All TacVar tick conversions (cntvct_el0 on aarch64, TSC on x86_64) go
 * through one runtime variable so a measured value (e.g. from
 * probe_counter_freq --nspt-out) can replace the ISA nominal rate:
 *
 *   g_tacvar_ns_per_tick
 *     - aarch64: loaded by tacvar_rate_init() from TACVAR_NSPT_FILE, else
 *       from the TACVAR_NSTP fallback;
 *     - x86_64:  set by tacvar_tsc_calibrate() (CPUID 0x15/0x16 or a
 *       CLOCK_MONOTONIC window) and mirrored here, so the same variable
 *       always holds the effective rate.
 *
 * Call tacvar_rate_init() from the tick-timer init path before any timed
 * region. Missing file / TACVAR_NSPT_FILE unset is not an error: the
 * nominal fallback is used.
 */
#ifndef TACVAR_TIMER_RATE_H
#define TACVAR_TIMER_RATE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Effective ns per hardware-counter tick (0 until initialized). */
extern double g_tacvar_ns_per_tick;

/**
 * @brief Load the measured ns-per-tick rate.
 * @param nspt_file  File written by probe_counter_freq --nspt-out (may be
 *                   NULL; the TACVAR_NSPT_FILE macro is tried next).
 * @param fallback   Nominal value used when no file is available
 *                   (e.g. TACVAR_NSTP on aarch64).
 * @return 0 on success (fallback or file), -1 only if the file exists but
 *         cannot be parsed (caller may then decide to abort).
 */
int tacvar_rate_init(const char *nspt_file, double fallback);

/**
 * @brief Parse a nspt file (first non-comment line) into the variable.
 * @return 0 on success, -1 if the value is missing or non-positive.
 */
int tacvar_rate_load_file(const char *nspt_file);

#ifdef __cplusplus
}
#endif

#endif /* TACVAR_TIMER_RATE_H */