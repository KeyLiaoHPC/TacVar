/**
 * @file counter_papi.h
 * @brief PAPI_read counter backend.
 */
#ifndef TACVAR_COUNTER_PAPI_H
#define TACVAR_COUNTER_PAPI_H

#include <stdint.h>

int  tacvar_counter_papi_init(const char *const *names, int n);
void tacvar_counter_papi_fini(void);
void tacvar_counter_papi_read(uint64_t *values);

#ifndef TACVAR_COUNTER_COUNT
#error "TACVAR_COUNTER_COUNT required"
#endif

#if TACVAR_COUNTER_COUNT > 0
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    int _i; \
    for (_i = 0; _i < TACVAR_COUNTER_COUNT; _i++) \
        (dst)[_i] = (uint64_t)((long long)(stop)[_i] - (long long)(start)[_i]); \
} while (0)
#else
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { (void)(dst); (void)(start); (void)(stop); } while (0)
#endif

#endif
