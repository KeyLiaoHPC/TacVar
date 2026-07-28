/**
 * @file tacvar_npb.c
 * @brief NPB-MPI TacVar timers: total filtering, per-step intervals, MPI init.
 *
 * Official NPB timer_read() still accumulates every timer slot. Only
 * whole-kernel totals (and optional TACVAR_REGION_STEP intervals) are
 * appended to the TacVar event buffer.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_npb.h"
#include "tacvar_npb_regions.h"
#include "tacvar_measure.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdint.h>
#include <sched.h>
#include <unistd.h>
#include <libgen.h>
#include <errno.h>
#include <limits.h>

#ifdef TACVAR_HAS_GENERATED_CONFIG
#include "tacvar_generated_config.h"
#endif

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif
#ifndef TACVAR_ENABLE_PER_STEP_TIMING
#define TACVAR_ENABLE_PER_STEP_TIMING 0
#endif
#ifndef TACVAR_REGION_STEP
#define TACVAR_REGION_STEP 1000
#endif

static double elapsed[64];
static uint64_t timer_start_raw[64];
static int cpu_start_slot[64];
#if TACVAR_COUNTER_COUNT > 0
static uint64_t counter_start_slot[64][TACVAR_COUNTER_COUNT];
#endif

typedef struct tacvar_interval {
    uint64_t raw_start;
#if TACVAR_COUNTER_COUNT > 0
    uint64_t counter_start[TACVAR_COUNTER_COUNT];
#endif
    int cpu_start;
    int active;
} tacvar_interval_t;

static tacvar_interval_t step_interval;
static tacvar_interval_t kernel_interval;

static char g_bench[64];
static char g_class[8];
static char g_test_tag[64];
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

static int is_tacvar_total_region(int slot)
{
    if (!strcasecmp(g_bench, "is"))
        return slot == 0;
    if (!strcasecmp(g_bench, "dt"))
        return 0; /* DT uses dedicated all-rank kernel interval */
    return slot == 1;
}

static void interval_begin(tacvar_interval_t *iv)
{
    iv->cpu_start = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
    if (tacvar_is_ready())
        TACVAR_COUNTER_READ(iv->counter_start);
#endif
    TACVAR_TIMER_BEGIN(iv->raw_start);
    iv->active = 1;
}

static void interval_end_and_record(tacvar_interval_t *iv, int region_id,
                                    const char *source)
{
    uint64_t raw_stop = 0;
    int64_t elapsed_ns;
    int cpu_stop;
#if TACVAR_COUNTER_COUNT > 0
    uint64_t counter_stop[TACVAR_COUNTER_COUNT];
    uint64_t counter_delta[TACVAR_COUNTER_COUNT];
#endif

    if (!iv->active)
        return;
    TACVAR_TIMER_END(raw_stop);
#if TACVAR_COUNTER_COUNT > 0
    if (tacvar_is_ready()) {
        TACVAR_COUNTER_READ(counter_stop);
        TACVAR_COUNTER_DELTAS(counter_delta, iv->counter_start, counter_stop);
    }
#endif
    cpu_stop = sched_getcpu();
    elapsed_ns = TACVAR_TIMER_DELTA_NS(iv->raw_start, raw_stop);
    iv->active = 0;

    if (!tacvar_is_ready())
        return;
    tacvar_csv_write_simple(region_id, iv->raw_start, raw_stop, elapsed_ns,
                            iv->cpu_start, cpu_stop,
#if TACVAR_COUNTER_COUNT > 0
                            iv->counter_start, counter_stop, counter_delta
#else
                            NULL, NULL, NULL
#endif
                            , source);
}

void tacvar_npb_set_test_tag(const char *tag)
{
    if (!tag) {
        g_test_tag[0] = '\0';
        return;
    }
    snprintf(g_test_tag, sizeof(g_test_tag), "%s", tag);
}

int tacvar_npb_ensure_init(void)
{
    tacvar_context_t ctx;
    int rank = 0, nprocs = 1;
    char dirbuf[512];
    int len = 0;
    int init_rc = 0;

    if (tacvar_is_ready())
        return 0;

    set_identity();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    memset(&ctx, 0, sizeof(ctx));
    dirbuf[0] = '\0';
    ctx.suite = "npb-mpi";
    ctx.benchmark = g_bench;
    ctx.klass = g_class;
    ctx.test_tag = g_test_tag[0] ? g_test_tag : NULL;
    ctx.rank = rank;
    ctx.thread = 0;
    ctx.nprocs = nprocs;

    /* Rank 0 must not return before Bcasts — peers would hang. */
    if (rank == 0) {
        init_rc = tacvar_init(&ctx);
        if (init_rc == 0) {
            snprintf(dirbuf, sizeof(dirbuf), "%s",
                     tacvar_data_dir() ? tacvar_data_dir() : "");
            len = (int)strlen(dirbuf) + 1;
        }
    }
    MPI_Bcast(&init_rc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (len > 0)
        MPI_Bcast(dirbuf, len, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (init_rc != 0)
        return init_rc;
    if (rank != 0) {
        setenv("TACVAR_DATA_DIR", dirbuf, 1);
        if (tacvar_init(&ctx) != 0)
            return -1;
    }
    /* Always re-register: tacvar_init zeroes region metadata on re-init/fork. */
    tacvar_npb_register_regions(g_bench);
    return 0;
}

void tacvar_npb_prepare(uint64_t expected_steps)
{
    size_t required;
    int rank = 0;
    int rc;

    if (tacvar_npb_ensure_init() != 0 && !g_init_warned) {
        fprintf(stderr, "tacvar: NPB-MPI init failed; counters/CSV disabled\n");
        g_init_warned = 1;
    }
    if (!tacvar_is_ready())
        return;

    required = 1;
#if TACVAR_ENABLE_PER_STEP_TIMING == 1
    if (expected_steps > SIZE_MAX - required) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        fprintf(stderr,
                "tacvar: rank %d event count overflow (steps=%llu)\n",
                rank, (unsigned long long)expected_steps);
        MPI_Abort(MPI_COMM_WORLD, ENOMEM);
    }
    required += (size_t)expected_steps;
#else
    (void)expected_steps;
#endif
    rc = tacvar_events_reserve(required);
    if (rc != 0) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        fprintf(stderr,
                "tacvar: rank %d cannot reserve %zu events (~%zu bytes): %s\n",
                rank, required, required * sizeof(uint64_t) * 8,
                strerror(-rc));
        MPI_Abort(MPI_COMM_WORLD, ENOMEM);
    }
}

void tacvar_npb_timer_clear(int n)
{
    if (n < 0 || n >= 64)
        return;
    elapsed[n] = 0.0;
}

void tacvar_npb_timer_start(int n)
{
    if (n < 0 || n >= 64)
        return;
    if (!tacvar_is_ready()) {
        if (tacvar_npb_ensure_init() != 0 && !g_init_warned) {
            fprintf(stderr, "tacvar: NPB-MPI init failed; counters/CSV disabled\n");
            g_init_warned = 1;
        }
    }
    cpu_start_slot[n] = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
    if (tacvar_is_ready() && is_tacvar_total_region(n))
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

    if (n < 0 || n >= 64)
        return;

    TACVAR_TIMER_END(timer_stop_raw);
    cpu_stop = sched_getcpu();
    elapsed_ns = TACVAR_TIMER_DELTA_NS(timer_start_raw[n], timer_stop_raw);
    elapsed[n] += (double)elapsed_ns * 1e-9;

    /* Always update NPB elapsed[]; only totals go to TacVar CSV. */
    if (!tacvar_is_ready() || !is_tacvar_total_region(n))
        return;

#if TACVAR_COUNTER_COUNT > 0
    TACVAR_COUNTER_READ(counter_stop);
    TACVAR_COUNTER_DELTAS(counter_delta, counter_start_slot[n], counter_stop);
#endif
    tacvar_csv_write_simple(n, timer_start_raw[n], timer_stop_raw, elapsed_ns,
                            cpu_start_slot[n], cpu_stop,
#if TACVAR_COUNTER_COUNT > 0
                            counter_start_slot[n], counter_stop, counter_delta
#else
                            NULL, NULL, NULL
#endif
                            , "npb-total");
}

double tacvar_npb_timer_read(int n)
{
    if (n < 0 || n >= 64)
        return 0.0;
    return elapsed[n];
}

void tacvar_npb_step_start(void)
{
#if TACVAR_ENABLE_PER_STEP_TIMING == 1
    if (!tacvar_is_ready()) {
        if (tacvar_npb_ensure_init() != 0)
            return;
    }
    interval_begin(&step_interval);
#else
    return;
#endif
}

void tacvar_npb_step_stop(void)
{
#if TACVAR_ENABLE_PER_STEP_TIMING == 1
    interval_end_and_record(&step_interval, TACVAR_REGION_STEP, "npb-step");
#else
    return;
#endif
}

void tacvar_npb_kernel_start(void)
{
    if (!tacvar_is_ready()) {
        if (tacvar_npb_ensure_init() != 0)
            return;
    }
    interval_begin(&kernel_interval);
}

void tacvar_npb_kernel_stop(void)
{
    interval_end_and_record(&kernel_interval, 0, "npb-kernel");
}
