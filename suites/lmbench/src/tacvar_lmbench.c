/**
 * @file tacvar_lmbench.c
 * @brief lmbench TacVar adapter: default global state + optional timeval sidecar.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_lmbench.h"
#include "tacvar_measure.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <libgen.h>
#include <sched.h>
#include <limits.h>

#ifdef TACVAR_HAS_GENERATED_CONFIG
#include "tacvar_generated_config.h"
#endif

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif

#define TACVAR_LMB_SIDECAR_MAX 8

typedef struct {
    struct timeval *key;
    uint64_t timer_start_raw;
    int cpu_start;
#if TACVAR_COUNTER_COUNT > 0
    uint64_t counter_start[TACVAR_COUNTER_COUNT];
#endif
    int used;
} tacvar_lmb_slot_t;

static uint64_t g_def_timer_start;
static int g_def_cpu_start;
#if TACVAR_COUNTER_COUNT > 0
static uint64_t g_def_counter_start[TACVAR_COUNTER_COUNT];
#endif
static tacvar_lmb_slot_t g_sidecar[TACVAR_LMB_SIDECAR_MAX];
static char g_bench[128];
static char g_tag[128];

static void encode_ns_to_tv(struct timeval *tv, uint64_t ns)
{
    if (!tv)
        return;
    tv->tv_sec = (time_t)(ns / 1000000000ULL);
    tv->tv_usec = (suseconds_t)((ns % 1000000000ULL) / 1000ULL);
}

static void set_identity(void)
{
    char path[PATH_MAX], base[256];
    ssize_t n;
    snprintf(g_bench, sizeof(g_bench), "lmbench");
    snprintf(g_tag, sizeof(g_tag), "%s",
             getenv("TACVAR_TEST_TAG") ? getenv("TACVAR_TEST_TAG") : "");
    n = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (n > 0) {
        path[n] = '\0';
        snprintf(base, sizeof(base), "%s", basename(path));
        snprintf(g_bench, sizeof(g_bench), "%s", base);
    }
}

void tacvar_lmbench_ensure_init(void)
{
    tacvar_context_t ctx;
    if (tacvar_is_ready())
        return;
    set_identity();
    memset(&ctx, 0, sizeof(ctx));
    ctx.suite = "lmbench";
    ctx.benchmark = g_bench;
    ctx.klass = "X";
    ctx.test_tag = g_tag;
    ctx.rank = 0;
    ctx.thread = 0;
    ctx.nprocs = 1;
    if (tacvar_init(&ctx) == 0) {
        int region_ids[1] = { 0 };
        int nlocs[1] = { 1 };
        const char *names[1] = { g_bench };
        (void)tacvar_write_timer_info(region_ids, nlocs, names, 1);
    }
}

static tacvar_lmb_slot_t *find_slot(struct timeval *key, int create)
{
    int i, free_i = -1;
    for (i = 0; i < TACVAR_LMB_SIDECAR_MAX; i++) {
        if (g_sidecar[i].used && g_sidecar[i].key == key)
            return &g_sidecar[i];
        if (!g_sidecar[i].used && free_i < 0)
            free_i = i;
    }
    if (!create || free_i < 0)
        return NULL;
    g_sidecar[free_i].used = 1;
    g_sidecar[free_i].key = key;
    return &g_sidecar[free_i];
}

void tacvar_lmbench_start(struct timeval *tv)
{
    tacvar_lmbench_ensure_init();
    if (tv == NULL) {
        g_def_cpu_start = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
        TACVAR_COUNTER_READ(g_def_counter_start);
#endif
        TACVAR_TIMER_BEGIN(g_def_timer_start);
        return;
    } else {
        tacvar_lmb_slot_t *s = find_slot(tv, 1);
        if (!s) {
            g_def_cpu_start = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
            TACVAR_COUNTER_READ(g_def_counter_start);
#endif
            TACVAR_TIMER_BEGIN(g_def_timer_start);
            encode_ns_to_tv(tv, g_def_timer_start);
            return;
        }
        s->cpu_start = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
        TACVAR_COUNTER_READ(s->counter_start);
#endif
        TACVAR_TIMER_BEGIN(s->timer_start_raw);
        encode_ns_to_tv(tv, s->timer_start_raw);
    }
}

uint64_t tacvar_lmbench_stop(struct timeval *begin, struct timeval *end)
{
    uint64_t timer_stop_raw = 0;
    uint64_t timer_start_raw;
    int cpu_start, cpu_stop;
    int64_t elapsed_ns;
#if TACVAR_COUNTER_COUNT > 0
    uint64_t counter_stop[TACVAR_COUNTER_COUNT];
    uint64_t counter_delta[TACVAR_COUNTER_COUNT];
    uint64_t *counter_start_ptr;
#endif
    tacvar_lmb_slot_t *s = NULL;

    tacvar_lmbench_ensure_init();

    if (begin != NULL)
        s = find_slot(begin, 0);

    TACVAR_TIMER_END(timer_stop_raw);
#if TACVAR_COUNTER_COUNT > 0
    TACVAR_COUNTER_READ(counter_stop);
#endif
    cpu_stop = sched_getcpu();

    if (s) {
        timer_start_raw = s->timer_start_raw;
        cpu_start = s->cpu_start;
#if TACVAR_COUNTER_COUNT > 0
        counter_start_ptr = s->counter_start;
        TACVAR_COUNTER_DELTAS(counter_delta, s->counter_start, counter_stop);
#endif
        s->used = 0;
    } else {
        timer_start_raw = g_def_timer_start;
        cpu_start = g_def_cpu_start;
#if TACVAR_COUNTER_COUNT > 0
        counter_start_ptr = g_def_counter_start;
        TACVAR_COUNTER_DELTAS(counter_delta, g_def_counter_start, counter_stop);
#endif
    }

    elapsed_ns = TACVAR_TIMER_DELTA_NS(timer_start_raw, timer_stop_raw);
    if (end)
        encode_ns_to_tv(end, timer_stop_raw);
    if (begin)
        encode_ns_to_tv(begin, timer_start_raw);

    tacvar_csv_write_simple(0, 0, timer_start_raw, timer_stop_raw, elapsed_ns,
                            cpu_start, cpu_stop,
#if TACVAR_COUNTER_COUNT > 0
                            counter_start_ptr, counter_stop, counter_delta
#else
                            NULL, NULL, NULL
#endif
                            );

    if (elapsed_ns < 0)
        elapsed_ns = 0;
    return (uint64_t)(elapsed_ns / 1000);
}

uint64_t tacvar_lmbench_now_us(void)
{
    uint64_t raw = 0;
    tacvar_lmbench_ensure_init();
    TACVAR_TIMER_BEGIN(raw);
    /* For ns backends, convert; for TSC raw convert via delta from 0 approx */
    return raw / 1000ULL;
}
