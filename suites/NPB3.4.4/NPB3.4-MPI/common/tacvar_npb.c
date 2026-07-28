/**
 * @file tacvar_npb.c
 * @brief NPB-MPI TacVar timers + MPI rank/dir context.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_npb.h"
#include "tacvar_measure.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sched.h>
#include <unistd.h>
#include <libgen.h>

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
        /* Local failure only; no second collective (would hang peers). */
        if (tacvar_init(&ctx) != 0)
            return -1;
    }
    return 0;
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
            fprintf(stderr, "tacvar: NPB-MPI init failed; counters/CSV disabled\n");
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

    if (tacvar_is_ready()) {
        tacvar_csv_write_simple(n, timer_start_raw[n], timer_stop_raw, elapsed_ns,
                                cpu_start_slot[n], cpu_stop,
#if TACVAR_COUNTER_COUNT > 0
                                counter_start_slot[n], counter_stop, counter_delta
#else
                                NULL, NULL, NULL
#endif
                                );
    }
}

double tacvar_npb_timer_read(int n)
{
    return elapsed[n];
}
