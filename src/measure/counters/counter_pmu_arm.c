/**
 * @file counter_pmu_arm.c
 * @brief Configure ARMv8 PMU event types at EL0 (kmod must enable access).
 */
#include "counter_pmu_arm.h"
#include "../events/events.h"
#include "../include/tacvar_measure.h"
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <errno.h>

unsigned char g_tacvar_arm_is_cycle[12];
static int g_n;
static uint32_t g_enable_mask;

static void tacvar_arm_set_pmevtyper(unsigned idx, uint32_t code)
{
    switch (idx) {
    case 0:  __asm__ volatile("msr pmevtyper0_el0, %0" :: "r"((uint64_t)code)); break;
    case 1:  __asm__ volatile("msr pmevtyper1_el0, %0" :: "r"((uint64_t)code)); break;
    case 2:  __asm__ volatile("msr pmevtyper2_el0, %0" :: "r"((uint64_t)code)); break;
    case 3:  __asm__ volatile("msr pmevtyper3_el0, %0" :: "r"((uint64_t)code)); break;
    case 4:  __asm__ volatile("msr pmevtyper4_el0, %0" :: "r"((uint64_t)code)); break;
    case 5:  __asm__ volatile("msr pmevtyper5_el0, %0" :: "r"((uint64_t)code)); break;
    case 6:  __asm__ volatile("msr pmevtyper6_el0, %0" :: "r"((uint64_t)code)); break;
    case 7:  __asm__ volatile("msr pmevtyper7_el0, %0" :: "r"((uint64_t)code)); break;
    case 8:  __asm__ volatile("msr pmevtyper8_el0, %0" :: "r"((uint64_t)code)); break;
    case 9:  __asm__ volatile("msr pmevtyper9_el0, %0" :: "r"((uint64_t)code)); break;
    case 10: __asm__ volatile("msr pmevtyper10_el0, %0" :: "r"((uint64_t)code)); break;
    case 11: __asm__ volatile("msr pmevtyper11_el0, %0" :: "r"((uint64_t)code)); break;
    default: break;
    }
}

int tacvar_counter_pmu_arm_init(const char *const *names, int n)
{
    int i, gp = 0;
    uint64_t pmcr, en;
    uint32_t ovf = 0xFFFFFFFFu;

    g_n = n;
    g_enable_mask = 0;
    memset(g_tacvar_arm_is_cycle, 0, sizeof(g_tacvar_arm_is_cycle));

    for (i = 0; i < n; i++) {
        uint32_t code = tacvar_armv8_parse_event(names[i]);
        if (code == 0xFFFFFFFFu) {
            fprintf(stderr, "tacvar: unknown ARM event '%s'\n", names[i]);
            return -EINVAL;
        }
        if (code == 0x11 || !strcasecmp(names[i], "CPU_CYCLES") ||
            !strcasecmp(names[i], "CYCLES")) {
            g_tacvar_arm_is_cycle[i] = 1;
            g_enable_mask |= (1u << 31);
        } else {
            if (gp > 11) {
                fprintf(stderr, "tacvar: too many GP ARM counters\n");
                return -EINVAL;
            }
            tacvar_arm_set_pmevtyper((unsigned)gp, code);
            g_enable_mask |= (1u << gp);
            gp++;
        }
    }

    /* Enable PMU, reset counters, clear overflow, enable selected. */
    __asm__ volatile("mrs %0, pmcr_el0" : "=r"(pmcr));
    pmcr |= 0x7; /* E | P | C */
    __asm__ volatile("msr pmcr_el0, %0" :: "r"(pmcr) : "memory");
    __asm__ volatile("msr pmovsclr_el0, %0" :: "r"((uint64_t)ovf) : "memory");
    en = (uint64_t)g_enable_mask;
    __asm__ volatile("msr pmcntenset_el0, %0" :: "r"(en) : "memory");
    return 0;
}

void tacvar_counter_pmu_arm_fini(void)
{
    uint64_t clr = (uint64_t)g_enable_mask;
    if (clr)
        __asm__ volatile("msr pmcntenclr_el0, %0" :: "r"(clr) : "memory");
    g_enable_mask = 0;
    g_n = 0;
}
