/**
 * @file timer_papi_real_nsec.h
 * @brief PAPI_get_real_nsec timer; raw is nanoseconds.
 */
#ifndef TACVAR_TIMER_PAPI_REAL_NSEC_H
#define TACVAR_TIMER_PAPI_REAL_NSEC_H

#include <papi.h>
#include <stdint.h>

static inline int tacvar_timer_papi_real_nsec_init(void)
{
    int rc = PAPI_library_init(PAPI_VER_CURRENT);
    if (rc != PAPI_VER_CURRENT && rc > 0)
        return -1;
    if (rc < 0)
        return rc;
    return 0;
}

static inline void tacvar_timer_papi_real_nsec_fini(void)
{
    /* Leave library up if counters also use PAPI. */
}

static inline void tacvar_timer_papi_real_nsec_begin(uint64_t *raw)
{
    *raw = (uint64_t)PAPI_get_real_nsec();
}

static inline void tacvar_timer_papi_real_nsec_end(uint64_t *raw)
{
    *raw = (uint64_t)PAPI_get_real_nsec();
}

static inline int64_t tacvar_timer_papi_real_nsec_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)(e - b);
}

#endif
