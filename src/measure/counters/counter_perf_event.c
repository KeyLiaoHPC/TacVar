/**
 * @file counter_perf_event.c
 * @brief Linux perf_event_open group reader for TacVar.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "counter_perf_event.h"
#include "../include/tacvar_measure.h"
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <strings.h>

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif

#if TACVAR_COUNTER_COUNT > 0

static int g_fds[TACVAR_MAX_COUNTERS];
static int g_n;
static int g_group = -1;

static long sys_perf_event_open(struct perf_event_attr *attr,
                                pid_t pid, int cpu, int group_fd, unsigned long flags)
{
    return syscall(__NR_perf_event_open, attr, pid, cpu, group_fd, flags);
}

static int parse_hw_type(const char *name, __u64 *config)
{
    if (!strcasecmp(name, "cpu-cycles") || !strcasecmp(name, "cycles") ||
        !strcasecmp(name, "CPU_CYCLES")) {
        *config = PERF_COUNT_HW_CPU_CYCLES;
        return PERF_TYPE_HARDWARE;
    }
    if (!strcasecmp(name, "instructions") || !strcasecmp(name, "INSTRUCTIONS") ||
        !strcasecmp(name, "INST_RETIRED")) {
        *config = PERF_COUNT_HW_INSTRUCTIONS;
        return PERF_TYPE_HARDWARE;
    }
    if (!strcasecmp(name, "cache-references") || !strcasecmp(name, "CACHE_REFERENCES")) {
        *config = PERF_COUNT_HW_CACHE_REFERENCES;
        return PERF_TYPE_HARDWARE;
    }
    if (!strcasecmp(name, "cache-misses") || !strcasecmp(name, "CACHE_MISSES")) {
        *config = PERF_COUNT_HW_CACHE_MISSES;
        return PERF_TYPE_HARDWARE;
    }
    if (!strcasecmp(name, "branch-instructions") || !strcasecmp(name, "BRANCHES")) {
        *config = PERF_COUNT_HW_BRANCH_INSTRUCTIONS;
        return PERF_TYPE_HARDWARE;
    }
    if (!strcasecmp(name, "branch-misses") || !strcasecmp(name, "BR_MIS_PRED")) {
        *config = PERF_COUNT_HW_BRANCH_MISSES;
        return PERF_TYPE_HARDWARE;
    }
    if (!strcasecmp(name, "ref-cycles")) {
        *config = PERF_COUNT_HW_REF_CPU_CYCLES;
        return PERF_TYPE_HARDWARE;
    }
    return -1;
}

int tacvar_counter_perf_event_init(const char *const *names, int n)
{
    int i;
    g_n = n;
    g_group = -1;
    for (i = 0; i < TACVAR_MAX_COUNTERS; i++)
        g_fds[i] = -1;

    for (i = 0; i < n; i++) {
        struct perf_event_attr attr;
        __u64 cfg = 0;
        int type;
        int fd;
        memset(&attr, 0, sizeof(attr));
        attr.size = sizeof(attr);
        attr.disabled = 1;
        /* exclude_kernel required when perf_event_paranoid >= 2 (non-root). */
        attr.exclude_kernel = 1;
        attr.exclude_hv = 1;
        type = parse_hw_type(names[i], &cfg);
        if (type < 0) {
            /* raw: rNNN or 0xNNN as PERF_TYPE_RAW */
            char *end = NULL;
            unsigned long long v;
            if (names[i][0] == 'r' || names[i][0] == 'R')
                v = strtoull(names[i] + 1, &end, 16);
            else if (names[i][0] == '0' && (names[i][1] == 'x' || names[i][1] == 'X'))
                v = strtoull(names[i], &end, 16);
            else {
                fprintf(stderr, "tacvar: unknown perf event '%s'\n", names[i]);
                return -EINVAL;
            }
            if (!end || *end != '\0') {
                fprintf(stderr, "tacvar: bad perf raw event '%s'\n", names[i]);
                return -EINVAL;
            }
            type = PERF_TYPE_RAW;
            cfg = (__u64)v;
        }
        attr.type = (__u32)type;
        attr.config = cfg;
        fd = (int)sys_perf_event_open(&attr, 0, -1, g_group, 0);
        if (fd < 0) {
            fprintf(stderr, "tacvar: perf_event_open(%s) failed: %s\n",
                    names[i], strerror(errno));
            tacvar_counter_perf_event_fini();
            return -errno;
        }
        g_fds[i] = fd;
        if (g_group < 0)
            g_group = fd;
    }
    if (g_group >= 0) {
        ioctl(g_group, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
        ioctl(g_group, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    }
    return 0;
}

void tacvar_counter_perf_event_fini(void)
{
    int i;
    if (g_group >= 0)
        ioctl(g_group, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
    for (i = 0; i < TACVAR_MAX_COUNTERS; i++) {
        if (g_fds[i] >= 0) {
            close(g_fds[i]);
            g_fds[i] = -1;
        }
    }
    g_group = -1;
    g_n = 0;
}

void tacvar_counter_perf_event_read(uint64_t *values)
{
    int i;
    for (i = 0; i < g_n; i++) {
        uint64_t v = 0;
        if (read(g_fds[i], &v, sizeof(v)) != (ssize_t)sizeof(v))
            v = 0;
        values[i] = v;
    }
}

#else /* TACVAR_COUNTER_COUNT == 0 */

int tacvar_counter_perf_event_init(const char *const *names, int n)
{
    (void)names; (void)n; return 0;
}
void tacvar_counter_perf_event_fini(void) {}
void tacvar_counter_perf_event_read(uint64_t *values) { (void)values; }

#endif
