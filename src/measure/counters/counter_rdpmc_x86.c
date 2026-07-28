/**
 * @file counter_rdpmc_x86.c
 * @brief Configure x86 GP PMCs via ph_enable_pmu sysfs; hot path uses RDPMC.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "counter_rdpmc_x86.h"
#include "../events/events.h"
#include "../include/tacvar_measure.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#define PH_SYSFS_CFG   "/sys/module/ph_enable_pmu/config"
#define PH_SYSFS_CNT   "/sys/module/ph_enable_pmu/counts"
#define PH_SYSFS_PCE   "/sys/module/ph_enable_pmu/cr4.pce"
#define PH_FIXED_SLOTS 3

static int g_cfg_fd = -1;
static int g_cnt_fd = -1;
static int g_n;

int tacvar_counter_rdpmc_x86_init(const char *const *names, int n)
{
    char pce[8];
    int i, fd;
    uint64_t codes[TACVAR_MAX_COUNTERS];
    uint64_t zeros[TACVAR_MAX_COUNTERS];
    ssize_t nw = -1;

    g_n = n;
    memset(zeros, 0, sizeof(zeros));

    fd = open(PH_SYSFS_PCE, O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "tacvar: cannot open %s (is ph_enable_pmu loaded?)\n",
                PH_SYSFS_PCE);
        return -ENOENT;
    }
    if (read(fd, pce, sizeof(pce) - 1) <= 0 || pce[0] != '1') {
        close(fd);
        fprintf(stderr, "tacvar: CR4.PCE not enabled\n");
        return -EPERM;
    }
    close(fd);

    for (i = 0; i < n; i++) {
        codes[i] = tacvar_x86_parse_event(names[i]);
        if (codes[i] == 0) {
            fprintf(stderr, "tacvar: unknown x86 event '%s'\n", names[i]);
            return -EINVAL;
        }
    }

    g_cfg_fd = open(PH_SYSFS_CFG, O_RDWR);
    g_cnt_fd = open(PH_SYSFS_CNT, O_RDWR);
    if (g_cfg_fd < 0 || g_cnt_fd < 0) {
        fprintf(stderr, "tacvar: cannot open ph_enable_pmu sysfs\n");
        tacvar_counter_rdpmc_x86_fini();
        return -ENOENT;
    }

    /* GP counters start at sysfs index 3 (after 3 fixed). */
    if (lseek(g_cfg_fd, (off_t)(PH_FIXED_SLOTS * sizeof(uint64_t)), SEEK_SET) < 0 ||
        (nw = write(g_cfg_fd, codes, (size_t)n * sizeof(uint64_t))) !=
            (ssize_t)((size_t)n * sizeof(uint64_t))) {
        fprintf(stderr, "tacvar: write config failed (%zd)\n", nw);
        tacvar_counter_rdpmc_x86_fini();
        return -EIO;
    }
    if (lseek(g_cnt_fd, (off_t)(PH_FIXED_SLOTS * sizeof(uint64_t)), SEEK_SET) < 0 ||
        write(g_cnt_fd, zeros, (size_t)n * sizeof(uint64_t)) !=
            (ssize_t)((size_t)n * sizeof(uint64_t))) {
        fprintf(stderr, "tacvar: zero counts failed\n");
        tacvar_counter_rdpmc_x86_fini();
        return -EIO;
    }
    return 0;
}

void tacvar_counter_rdpmc_x86_fini(void)
{
    if (g_cfg_fd >= 0) { close(g_cfg_fd); g_cfg_fd = -1; }
    if (g_cnt_fd >= 0) { close(g_cnt_fd); g_cnt_fd = -1; }
    g_n = 0;
}
