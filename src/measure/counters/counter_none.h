/**
 * @file counter_none.h
 * @brief No-op counter backend (compile-time empty macros).
 */
#ifndef TACVAR_COUNTER_NONE_H
#define TACVAR_COUNTER_NONE_H

#include <stdint.h>

static inline int tacvar_counter_none_init(const char *const *names, int n)
{
    (void)names;
    (void)n;
    return 0;
}
static inline void tacvar_counter_none_fini(void) {}
static inline void tacvar_counter_none_read(uint64_t *values)
{
    (void)values;
}

#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { (void)(dst); (void)(start); (void)(stop); } while (0)

#endif
