/**
 * @file tacvar_npb.c
 * @brief NPB-OMP TacVar timers with OpenMP thread-private state.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_npb.h"
#include "tacvar_measure.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sched.h>
#include <unistd.h>
#include <libgen.h>
#include <dlfcn.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef TACVAR_HAS_GENERATED_CONFIG
#include "tacvar_generated_config.h"
#endif

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif

static double elapsed[64];
static uint64_t timer_start_raw[64];
static int cpu_start_slot[64];
#if TACVAR_COUNTER_COUNT > 0
static uint64_t counter_start_slot[64][TACVAR_COUNTER_COUNT];
#endif

#ifdef _OPENMP
#pragma omp threadprivate(elapsed, timer_start_raw, cpu_start_slot)
#if TACVAR_COUNTER_COUNT > 0
#pragma omp threadprivate(counter_start_slot)
#endif
#endif

/* Per-timer-ID source location cache (resolved via dladdr). */
static char timer_source[64][128];

/**
 * Resolve a return address to a human-readable symbol string.
 * Result is cached per timer ID so dladdr runs only once per timer slot.
 */
static const char *resolve_source(int n, void *retaddr)
{
    char *buf = timer_source[n];
    if (buf[0] != '\0')
        return buf;

    buf[0] = '\0';
    if (!retaddr)
        return buf;

    Dl_info info;
    if (dladdr(retaddr, &info) && info.dli_sname) {
        snprintf(buf, sizeof(timer_source[0]), "%s", info.dli_sname);
    }
    return buf;
}

static char g_bench[64];
static char g_class[8];
static int g_identity_set;
static int g_init_warned;

static void set_identity(void)
{
    char path[512], base[256], *b, *p;
    ssize_t n;
    if (g_identity_set)
        return;
    snprintf(g_bench, sizeof(g_bench), "unknown");
    snprintf(g_class, sizeof(g_class), "U");
    n = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (n > 0) {
        path[n] = '\0';
        snprintf(base, sizeof(base), "%s", basename(path));
        b = base;
        p = strchr(b, '.');
        if (p) {
            *p = '\0';
            snprintf(g_bench, sizeof(g_bench), "%s", b);
            if (p[1]) {
                g_class[0] = p[1];
                g_class[1] = '\0';
            }
        } else {
            snprintf(g_bench, sizeof(g_bench), "%s", base);
        }
    }
    g_identity_set = 1;
}

int tacvar_npb_ensure_init(void)
{
    tacvar_context_t ctx;
    int thread = 0;
    int rc = 0;

    if (tacvar_is_ready())
        return 0;

#ifdef _OPENMP
#pragma omp critical(tacvar_npb_init)
#endif
    {
        if (!tacvar_is_ready()) {
            set_identity();
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif
            memset(&ctx, 0, sizeof(ctx));
            ctx.suite = "npb-omp";
            ctx.benchmark = g_bench;
            ctx.klass = g_class;
            ctx.rank = 0;
            ctx.thread = thread;
            ctx.nprocs = 1;
#ifdef _OPENMP
            ctx.nprocs = omp_get_max_threads();
#endif
            rc = tacvar_init(&ctx);
        }
    }
    return tacvar_is_ready() ? 0 : (rc != 0 ? rc : -1);
}

void tacvar_npb_timer_clear(int n)
{
    elapsed[n] = 0.0;
}

void tacvar_npb_timer_start(int n)
{
    if (!tacvar_is_ready()) {
        if (tacvar_npb_ensure_init() != 0 && !g_init_warned) {
            /* native timers need no TIMER_INIT; TSC may be wrong if init never ran. */
            fprintf(stderr, "tacvar: NPB-OMP init failed; counters/CSV disabled\n");
            g_init_warned = 1;
        }
    }
    cpu_start_slot[n] = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
    if (tacvar_is_ready())
        TACVAR_COUNTER_READ(counter_start_slot[n]);
#endif
    TACVAR_TIMER_BEGIN(timer_start_raw[n]);
}

void tacvar_npb_timer_stop(int n)
{
    uint64_t timer_stop_raw = 0;
    int64_t elapsed_ns;
    int cpu_stop;
#if TACVAR_COUNTER_COUNT > 0
    uint64_t counter_stop[TACVAR_COUNTER_COUNT];
    uint64_t counter_delta[TACVAR_COUNTER_COUNT];
#endif
    const char *source = "";

    TACVAR_TIMER_END(timer_stop_raw);
#if TACVAR_COUNTER_COUNT > 0
    if (tacvar_is_ready()) {
        TACVAR_COUNTER_READ(counter_stop);
        TACVAR_COUNTER_DELTAS(counter_delta, counter_start_slot[n], counter_stop);
    }
#endif
    cpu_stop = sched_getcpu();
    elapsed_ns = TACVAR_TIMER_DELTA_NS(timer_start_raw[n], timer_stop_raw);
    elapsed[n] += (double)elapsed_ns * 1e-9;

    /* Resolve caller source location (cached per timer ID). */
    source = resolve_source(n, __builtin_return_address(0));

    if (tacvar_is_ready()) {
        tacvar_csv_write_simple(n, timer_start_raw[n], timer_stop_raw, elapsed_ns,
                                cpu_start_slot[n], cpu_stop,
#if TACVAR_COUNTER_COUNT > 0
                                counter_start_slot[n], counter_stop, counter_delta
#else
                                NULL, NULL, NULL
#endif
                                , source);
    }
}

double tacvar_npb_timer_read(int n)
{
    return elapsed[n];
}
