/**
 * @file counter_perf_event.h
 * @brief perf_event_open counter backend (hot-path read via function call).
 */
#ifndef TACVAR_COUNTER_PERF_EVENT_H
#define TACVAR_COUNTER_PERF_EVENT_H

#include <stdint.h>

int  tacvar_counter_perf_event_init(const char *const *names, int n);
void tacvar_counter_perf_event_fini(void);
void tacvar_counter_perf_event_read(uint64_t *values);

#ifndef TACVAR_COUNTER_COUNT
#error "TACVAR_COUNTER_COUNT required"
#endif

#if TACVAR_COUNTER_COUNT > 0
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    int _i; \
    for (_i = 0; _i < TACVAR_COUNTER_COUNT; _i++) \
        (dst)[_i] = (stop)[_i] - (start)[_i]; \
} while (0)
#else
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { (void)(dst); (void)(start); (void)(stop); } while (0)
#endif

#endif
