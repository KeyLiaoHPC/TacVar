/**
 * @file counter_rdpmc_x86.h
 * @brief Hot-path RDPMC reads; count is fixed at generate time.
 *
 * GP counters use RDPMC indices 0..N-1 (programmed into PERFEVTSEL via kmod
 * slots starting at index 3). Fixed counters use 0x40000000|idx if configured
 * by name FIXED0/FIXED1/FIXED2.
 */
#ifndef TACVAR_COUNTER_RDPMC_X86_H
#define TACVAR_COUNTER_RDPMC_X86_H

#include <stdint.h>

int  tacvar_counter_rdpmc_x86_init(const char *const *names, int n);
void tacvar_counter_rdpmc_x86_fini(void);

static inline uint64_t tacvar_rdpmc(unsigned idx)
{
    unsigned lo, hi;
    __asm__ volatile("rdpmc" : "=a"(lo), "=d"(hi) : "c"(idx) : "memory");
    return ((uint64_t)hi << 32) | (uint64_t)lo;
}

#ifndef TACVAR_COUNTER_COUNT
#error "TACVAR_COUNTER_COUNT required"
#endif

#if TACVAR_COUNTER_COUNT == 0
static inline void tacvar_counter_rdpmc_x86_read(uint64_t *values)
{
    (void)values;
}
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { (void)(dst); (void)(start); (void)(stop); } while (0)

#elif TACVAR_COUNTER_COUNT == 1
static inline void tacvar_counter_rdpmc_x86_read(uint64_t *values)
{
    values[0] = tacvar_rdpmc(0);
}
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { (dst)[0] = (stop)[0] - (start)[0]; } while (0)

#elif TACVAR_COUNTER_COUNT == 2
static inline void tacvar_counter_rdpmc_x86_read(uint64_t *values)
{
    values[0] = tacvar_rdpmc(0);
    values[1] = tacvar_rdpmc(1);
}
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    (dst)[0] = (stop)[0] - (start)[0]; \
    (dst)[1] = (stop)[1] - (start)[1]; \
} while (0)

#elif TACVAR_COUNTER_COUNT == 3
static inline void tacvar_counter_rdpmc_x86_read(uint64_t *values)
{
    values[0] = tacvar_rdpmc(0);
    values[1] = tacvar_rdpmc(1);
    values[2] = tacvar_rdpmc(2);
}
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    (dst)[0] = (stop)[0] - (start)[0]; \
    (dst)[1] = (stop)[1] - (start)[1]; \
    (dst)[2] = (stop)[2] - (start)[2]; \
} while (0)

#elif TACVAR_COUNTER_COUNT == 4
static inline void tacvar_counter_rdpmc_x86_read(uint64_t *values)
{
    values[0] = tacvar_rdpmc(0);
    values[1] = tacvar_rdpmc(1);
    values[2] = tacvar_rdpmc(2);
    values[3] = tacvar_rdpmc(3);
}
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    int _i; for (_i = 0; _i < 4; _i++) (dst)[_i] = (stop)[_i] - (start)[_i]; \
} while (0)

#else
static inline void tacvar_counter_rdpmc_x86_read(uint64_t *values)
{
    int i;
    for (i = 0; i < TACVAR_COUNTER_COUNT; i++)
        values[i] = tacvar_rdpmc((unsigned)i);
}
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    int _i; for (_i = 0; _i < TACVAR_COUNTER_COUNT; _i++) \
        (dst)[_i] = (stop)[_i] - (start)[_i]; \
} while (0)
#endif

#endif
