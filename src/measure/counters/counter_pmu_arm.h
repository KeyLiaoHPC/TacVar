/**
 * @file counter_pmu_arm.h
 * @brief ARMv8 PMEVCNTR/PMCCNTR hot-path reads (compile-time count).
 */
#ifndef TACVAR_COUNTER_PMU_ARM_H
#define TACVAR_COUNTER_PMU_ARM_H

#include <stdint.h>

int  tacvar_counter_pmu_arm_init(const char *const *names, int n);
void tacvar_counter_pmu_arm_fini(void);

/* Which slots use cycle counter (PMCCNTR) vs GP (PMEVCNTR). Filled at init. */
extern unsigned char g_tacvar_arm_is_cycle[12];

static inline uint64_t tacvar_arm_read_pmccntr(void)
{
    uint64_t v;
    __asm__ volatile("mrs %0, pmccntr_el0" : "=r"(v) :: "memory");
    return v;
}

static inline uint64_t tacvar_arm_read_pmevcntr(unsigned idx)
{
    uint64_t v = 0;
    switch (idx) {
    case 0:  __asm__ volatile("mrs %0, pmevcntr0_el0"  : "=r"(v)); break;
    case 1:  __asm__ volatile("mrs %0, pmevcntr1_el0"  : "=r"(v)); break;
    case 2:  __asm__ volatile("mrs %0, pmevcntr2_el0"  : "=r"(v)); break;
    case 3:  __asm__ volatile("mrs %0, pmevcntr3_el0"  : "=r"(v)); break;
    case 4:  __asm__ volatile("mrs %0, pmevcntr4_el0"  : "=r"(v)); break;
    case 5:  __asm__ volatile("mrs %0, pmevcntr5_el0"  : "=r"(v)); break;
    case 6:  __asm__ volatile("mrs %0, pmevcntr6_el0"  : "=r"(v)); break;
    case 7:  __asm__ volatile("mrs %0, pmevcntr7_el0"  : "=r"(v)); break;
    case 8:  __asm__ volatile("mrs %0, pmevcntr8_el0"  : "=r"(v)); break;
    case 9:  __asm__ volatile("mrs %0, pmevcntr9_el0"  : "=r"(v)); break;
    case 10: __asm__ volatile("mrs %0, pmevcntr10_el0" : "=r"(v)); break;
    case 11: __asm__ volatile("mrs %0, pmevcntr11_el0" : "=r"(v)); break;
    default: break;
    }
    return v;
}

#ifndef TACVAR_COUNTER_COUNT
#error "TACVAR_COUNTER_COUNT required"
#endif

static inline void tacvar_counter_pmu_arm_read(uint64_t *values)
{
#if TACVAR_COUNTER_COUNT == 0
    (void)values;
#else
    int i;
    int gp = 0;
    for (i = 0; i < TACVAR_COUNTER_COUNT; i++) {
        if (g_tacvar_arm_is_cycle[i])
            values[i] = tacvar_arm_read_pmccntr();
        else {
            values[i] = tacvar_arm_read_pmevcntr((unsigned)gp);
            gp++;
        }
    }
#endif
}

#if TACVAR_COUNTER_COUNT > 0
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { \
    int _i; for (_i = 0; _i < TACVAR_COUNTER_COUNT; _i++) \
        (dst)[_i] = (stop)[_i] - (start)[_i]; \
} while (0)
#else
#define TACVAR_COUNTER_DELTAS(dst, start, stop) do { (void)(dst); (void)(start); (void)(stop); } while (0)
#endif

#endif
