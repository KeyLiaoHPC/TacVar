/**
 * @file tacvar_npb.h
 * @brief NPB-MPI TacVar timer/adapter API (C and Fortran bind(C)).
 *
 * Existing NPB timer_start/stop continue to accumulate elapsed[] for the
 * official NPB report. Only whole-kernel totals are appended to TacVar
 * event buffers. Optional per-step intervals use TACVAR_REGION_STEP.
 */
#ifndef TACVAR_NPB_H
#define TACVAR_NPB_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void tacvar_npb_timer_clear(int n);
void tacvar_npb_timer_start(int n);
void tacvar_npb_timer_stop(int n);
double tacvar_npb_timer_read(int n);
int tacvar_npb_ensure_init(void);

/** Preallocate event buffer (1 total + optional steps). Aborts on OOM. */
void tacvar_npb_prepare(uint64_t expected_steps);

/** Per-step interval (no-op when TACVAR_ENABLE_PER_STEP_TIMING=0). */
void tacvar_npb_step_start(void);
void tacvar_npb_step_stop(void);

/** DT all-rank whole-kernel interval (independent of NPB slot 0). */
void tacvar_npb_kernel_start(void);
void tacvar_npb_kernel_stop(void);

/** Optional test_tag (e.g. DT graph BH/WH/SH); call before ensure_init. */
void tacvar_npb_set_test_tag(const char *tag);

#ifdef __cplusplus
}
#endif
#endif
