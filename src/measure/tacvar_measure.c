/**
 * @file tacvar_measure.c
 * @brief Non-hot-path init/fini and fork-safe measurement state.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#ifdef TACVAR_HAS_GENERATED_CONFIG
#include "tacvar_generated_config.h"
#endif

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif
#ifndef TACVAR_TIMER_NAME
#define TACVAR_TIMER_NAME "unknown"
#endif
#ifndef TACVAR_COUNTER_BACKEND_NAME
#define TACVAR_COUNTER_BACKEND_NAME "none"
#endif
#ifndef TACVAR_OUTPUT_ROOT_DEFAULT
#define TACVAR_OUTPUT_ROOT_DEFAULT "."
#endif

tacvar_state_t g_tacvar;

uint64_t tacvar_modular_delta(uint64_t start, uint64_t stop, unsigned width_bits)
{
    if (width_bits == 0 || width_bits >= 64)
        return stop - start;
    {
        uint64_t mask = (width_bits == 64) ? ~0ULL : ((1ULL << width_bits) - 1ULL);
        return (stop - start) & mask;
    }
}

static void copy_context(tacvar_context_t *dst, const tacvar_context_t *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
        return;
    *dst = *src;
}

static void fill_names_from_generated(tacvar_state_t *st)
{
    int i;
    snprintf(st->timer_name, sizeof(st->timer_name), "%s", TACVAR_TIMER_NAME);
    snprintf(st->counter_backend, sizeof(st->counter_backend), "%s",
             TACVAR_COUNTER_BACKEND_NAME);
    st->n_counters = TACVAR_COUNTER_COUNT;
#if TACVAR_COUNTER_COUNT > 0
    for (i = 0; i < TACVAR_COUNTER_COUNT; i++) {
        snprintf(st->counter_names[i], sizeof(st->counter_names[i]), "%s",
                 TACVAR_COUNTER_NAME_LIST[i]);
    }
#else
    (void)i;
#endif
}

int tacvar_init(const tacvar_context_t *ctx)
{
    pid_t pid = getpid();
    int rc;

    if (g_tacvar.initialized && g_tacvar.owner_pid == pid)
        return 0;

    if (g_tacvar.initialized && g_tacvar.owner_pid != pid) {
        /* Inherited state after fork — discard without flushing parent. */
        g_tacvar.csv_fp = NULL;
        g_tacvar.events = NULL;
        g_tacvar.event_count = 0;
        g_tacvar.event_capacity = 0;
        g_tacvar.initialized = 0;
#ifdef TACVAR_HAS_GENERATED_CONFIG
        TACVAR_COUNTER_FINI();
#endif
    }

    /* Leftover buffer from a prior failed flush in this process. */
    if (!g_tacvar.initialized && g_tacvar.events != NULL) {
        free(g_tacvar.events);
        g_tacvar.events = NULL;
        g_tacvar.event_count = 0;
        g_tacvar.event_capacity = 0;
    }

    memset(&g_tacvar, 0, sizeof(g_tacvar));
    copy_context(&g_tacvar.ctx, ctx);
    fill_names_from_generated(&g_tacvar);

    rc = tacvar_prepare_data_dir(&g_tacvar, TACVAR_OUTPUT_ROOT_DEFAULT);
    if (rc != 0) {
        fprintf(stderr, "tacvar: prepare data dir failed (%d)\n", rc);
        return rc;
    }

#ifdef TACVAR_HAS_GENERATED_CONFIG
    rc = TACVAR_TIMER_INIT();
    if (rc != 0) {
        fprintf(stderr, "tacvar: timer init failed (%d)\n", rc);
        return rc;
    }
    rc = TACVAR_COUNTER_INIT();
    if (rc != 0) {
        fprintf(stderr, "tacvar: counter init failed (%d)\n", rc);
        TACVAR_TIMER_FINI();
        return rc;
    }
#endif

    rc = tacvar_csv_open(&g_tacvar);
    if (rc != 0) {
        fprintf(stderr, "tacvar: csv open failed (%d)\n", rc);
#ifdef TACVAR_HAS_GENERATED_CONFIG
        TACVAR_COUNTER_FINI();
        TACVAR_TIMER_FINI();
#endif
        return rc;
    }

    g_tacvar.owner_pid = pid;
    g_tacvar.initialized = 1;
    atexit(tacvar_fini);
    return 0;
}

void tacvar_fini(void)
{
    if (!g_tacvar.initialized)
        return;
    if (g_tacvar.owner_pid != getpid())
        return;
    (void)tacvar_region_info_write();
    (void)tacvar_csv_close(&g_tacvar);
#ifdef TACVAR_HAS_GENERATED_CONFIG
    TACVAR_COUNTER_FINI();
    TACVAR_TIMER_FINI();
#endif
    g_tacvar.initialized = 0;
}

const char *tacvar_data_dir(void)
{
    return g_tacvar.initialized ? g_tacvar.data_dir : NULL;
}

int tacvar_is_ready(void)
{
    return g_tacvar.initialized && g_tacvar.owner_pid == getpid();
}
