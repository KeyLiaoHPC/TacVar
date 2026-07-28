/**
 * @file timer_omp_get_wtime.h
 * @brief OpenMP wall timer with original NPB-OMP wtime fallback.
 */
#ifndef TACVAR_TIMER_OMP_GET_WTIME_H
#define TACVAR_TIMER_OMP_GET_WTIME_H

#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#else
void wtime(double *t);
#endif

static inline int tacvar_timer_omp_get_wtime_init(void) { return 0; }
static inline void tacvar_timer_omp_get_wtime_fini(void) {}

static inline void tacvar_timer_omp_get_wtime_begin(uint64_t *raw)
{
#ifdef _OPENMP
    *raw = (uint64_t)(omp_get_wtime() * 1e9);
#else
    double t;
    wtime(&t);
    *raw = (uint64_t)(t * 1e9);
#endif
}

static inline void tacvar_timer_omp_get_wtime_end(uint64_t *raw)
{
#ifdef _OPENMP
    *raw = (uint64_t)(omp_get_wtime() * 1e9);
#else
    double t;
    wtime(&t);
    *raw = (uint64_t)(t * 1e9);
#endif
}

static inline int64_t tacvar_timer_omp_get_wtime_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)(e - b);
}

#endif
