/**
 * @file timer_util.c
 * @brief TSC frequency calibration for x86 TacVar timers.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "timer_util.h"
#include <stdio.h>
#include <unistd.h>

#ifdef __x86_64__

double g_tacvar_ns_per_tsc = 0.0;

static int tacvar_cpuid_tsc_hz(uint64_t *hz_out)
{
    unsigned eax, ebx, ecx, edx;
    /* Leaf 0x15: TSC/core crystal ratio */
    eax = 0x15;
    ecx = 0;
    __asm__ volatile("cpuid"
                     : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                     : "a"(eax), "c"(ecx));
    if (eax != 0 && ebx != 0 && ecx != 0) {
        *hz_out = ((uint64_t)ecx * (uint64_t)ebx) / (uint64_t)eax;
        return 0;
    }
    /* Leaf 0x16: base frequency in MHz (eax) */
    eax = 0x16;
    ecx = 0;
    __asm__ volatile("cpuid"
                     : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                     : "a"(eax), "c"(ecx));
    if (eax != 0) {
        *hz_out = (uint64_t)eax * 1000000ULL;
        return 0;
    }
    return -1;
}

static int tacvar_calibrate_with_cgt(uint64_t *hz_out)
{
    struct timespec t0, t1;
    unsigned hi0, lo0, hi1, lo1;
    uint64_t c0, c1;
    double secs;

    clock_gettime(CLOCK_MONOTONIC, &t0);
    __asm__ volatile("rdtsc" : "=a"(lo0), "=d"(hi0) :: "memory");
    usleep(200000); /* 200 ms */
    __asm__ volatile("rdtsc" : "=a"(lo1), "=d"(hi1) :: "memory");
    clock_gettime(CLOCK_MONOTONIC, &t1);

    c0 = tacvar_join_edx_eax(hi0, lo0);
    c1 = tacvar_join_edx_eax(hi1, lo1);
    secs = (double)(t1.tv_sec - t0.tv_sec)
         + (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
    if (secs <= 0.0 || c1 <= c0)
        return -1;
    *hz_out = (uint64_t)((double)(c1 - c0) / secs);
    return 0;
}

int tacvar_tsc_calibrate(void)
{
    uint64_t hz = 0;
    if (tacvar_cpuid_tsc_hz(&hz) != 0) {
        if (tacvar_calibrate_with_cgt(&hz) != 0 || hz == 0)
            return -1;
    }
    g_tacvar_ns_per_tsc = 1e9 / (double)hz;
    return 0;
}

#endif /* __x86_64__ */
