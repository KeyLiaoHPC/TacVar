/**
 * @file tacvar_npb.c
 * @brief NPB-OMP TacVar timers with OpenMP thread-private state + timer_info.csv.
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
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef TACVAR_HAS_GENERATED_CONFIG
#include "tacvar_generated_config.h"
#endif

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif

#define TACVAR_NPB_SLOTS 64

static double elapsed[TACVAR_NPB_SLOTS];
static uint64_t timer_start_raw[TACVAR_NPB_SLOTS];
static int cpu_start_slot[TACVAR_NPB_SLOTS];
static int max_loc[TACVAR_NPB_SLOTS];
#if TACVAR_COUNTER_COUNT > 0
static uint64_t counter_start_slot[TACVAR_NPB_SLOTS][TACVAR_COUNTER_COUNT];
#endif

#ifdef _OPENMP
#pragma omp threadprivate(elapsed, timer_start_raw, cpu_start_slot)
#if TACVAR_COUNTER_COUNT > 0
#pragma omp threadprivate(counter_start_slot)
#endif
#endif

static char g_bench[64];
static char g_class[8];
static int g_identity_set;
static int g_init_warned;
static int g_timer_info_registered;

typedef struct {
    int region_id;
    const char *name;
} tacvar_npb_name_t;

static const tacvar_npb_name_t names_is[] = {
    {0, "total"}, {1, "rank"}, {2, "rcomm"}, {3, "verify"},
};
static const tacvar_npb_name_t names_cg[] = {
    {1, "init"}, {2, "benchmk"}, {3, "conjgd"},
};
static const tacvar_npb_name_t names_mg[] = {
    {1, "init"}, {2, "benchmk"}, {3, "mg3P"}, {4, "psinv"}, {5, "resid"},
    {6, "resid2"}, {7, "rprj3"}, {8, "interp"}, {9, "norm2"}, {10, "comm3"},
};
static const tacvar_npb_name_t names_ep[] = {
    {1, "total"}, {2, "guv"}, {3, "randn"},
};
static const tacvar_npb_name_t names_bt[] = {
    {1, "total"}, {2, "rhsx"}, {3, "rhsy"}, {4, "rhsz"}, {5, "rhs"},
    {6, "xsolve"}, {7, "ysolve"}, {8, "zsolve"}, {9, "redist1"},
    {10, "solsubs"}, {11, "add"},
};
static const tacvar_npb_name_t names_sp[] = {
    {1, "total"}, {2, "rhsx"}, {3, "rhsy"}, {4, "rhsz"}, {5, "rhs"},
    {6, "xsolve"}, {7, "ysolve"}, {8, "zsolve"}, {9, "redist1"},
    {10, "redist2"}, {11, "txinvr"}, {12, "pinvr"}, {13, "ninvr"},
    {14, "tzetar"}, {15, "add"},
};
static const tacvar_npb_name_t names_lu[] = {
    {1, "total"}, {2, "rhsx"}, {3, "rhsy"}, {4, "rhsz"}, {5, "rhs"},
    {6, "jacld"}, {7, "blts"}, {8, "jacu"}, {9, "buts"}, {10, "add"},
    {11, "l2norm"},
};
static const tacvar_npb_name_t names_ft[] = {
    {1, "total"}, {2, "setup"}, {3, "fft"}, {4, "evolve"}, {5, "checksum"},
    {6, "fftx"}, {7, "ffty"}, {8, "fftz"},
};
static const tacvar_npb_name_t names_ua[] = {
    {1, "total"}, {2, "init"}, {3, "convect"}, {4, "transfb_c"},
    {5, "diffusion"}, {6, "transf"}, {7, "transfb"}, {8, "adaptation"},
    {9, "transf+b"}, {10, "add2"},
};
static const tacvar_npb_name_t names_dc[] = {
    {0, "total"},
};

static const tacvar_npb_name_t *lookup_names(int *n_out)
{
    if (!strcmp(g_bench, "is")) { *n_out = (int)(sizeof(names_is)/sizeof(names_is[0])); return names_is; }
    if (!strcmp(g_bench, "cg")) { *n_out = (int)(sizeof(names_cg)/sizeof(names_cg[0])); return names_cg; }
    if (!strcmp(g_bench, "mg")) { *n_out = (int)(sizeof(names_mg)/sizeof(names_mg[0])); return names_mg; }
    if (!strcmp(g_bench, "ep")) { *n_out = (int)(sizeof(names_ep)/sizeof(names_ep[0])); return names_ep; }
    if (!strcmp(g_bench, "bt")) { *n_out = (int)(sizeof(names_bt)/sizeof(names_bt[0])); return names_bt; }
    if (!strcmp(g_bench, "sp")) { *n_out = (int)(sizeof(names_sp)/sizeof(names_sp[0])); return names_sp; }
    if (!strcmp(g_bench, "lu")) { *n_out = (int)(sizeof(names_lu)/sizeof(names_lu[0])); return names_lu; }
    if (!strcmp(g_bench, "ft")) { *n_out = (int)(sizeof(names_ft)/sizeof(names_ft[0])); return names_ft; }
    if (!strcmp(g_bench, "ua")) { *n_out = (int)(sizeof(names_ua)/sizeof(names_ua[0])); return names_ua; }
    if (!strcmp(g_bench, "dc")) { *n_out = (int)(sizeof(names_dc)/sizeof(names_dc[0])); return names_dc; }
    *n_out = 0;
    return NULL;
}

static void note_loc(int n, int loc_id)
{
    if (n < 0 || n >= TACVAR_NPB_SLOTS)
        return;
#ifdef _OPENMP
#pragma omp critical(tacvar_npb_maxloc)
#endif
    {
        if (loc_id > max_loc[n])
            max_loc[n] = loc_id;
    }
}

static void write_timer_info_atexit(void)
{
    const tacvar_npb_name_t *names;
    int nnames, i, n;
    int region_ids[TACVAR_NPB_SLOTS];
    int nlocs[TACVAR_NPB_SLOTS];
    const char *nameptrs[TACVAR_NPB_SLOTS];
    char fallback[TACVAR_NPB_SLOTS][16];
    int have_name[TACVAR_NPB_SLOTS];

    if (!tacvar_is_ready())
        return;

    memset(have_name, 0, sizeof(have_name));
    names = lookup_names(&nnames);
    for (i = 0; i < nnames; i++) {
        int rid = names[i].region_id;
        if (rid >= 0 && rid < TACVAR_NPB_SLOTS) {
            have_name[rid] = 1;
            if (max_loc[rid] < 1)
                max_loc[rid] = 1;
        }
    }

    n = 0;
    for (i = 0; i < TACVAR_NPB_SLOTS; i++) {
        if (max_loc[i] <= 0 && !have_name[i])
            continue;
        region_ids[n] = i;
        nlocs[n] = max_loc[i] > 0 ? max_loc[i] : 1;
        nameptrs[n] = NULL;
        if (names) {
            int j;
            for (j = 0; j < nnames; j++) {
                if (names[j].region_id == i) {
                    nameptrs[n] = names[j].name;
                    break;
                }
            }
        }
        if (!nameptrs[n]) {
            snprintf(fallback[n], sizeof(fallback[n]), "r%d", i);
            nameptrs[n] = fallback[n];
        }
        n++;
    }
    if (n > 0)
        (void)tacvar_write_timer_info(region_ids, nlocs, nameptrs, n);
}

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
            if (rc == 0 && !g_timer_info_registered) {
                atexit(write_timer_info_atexit);
                g_timer_info_registered = 1;
            }
        }
    }
    return tacvar_is_ready() ? 0 : (rc != 0 ? rc : -1);
}

void tacvar_npb_timer_clear(int n)
{
    if (n < 0 || n >= TACVAR_NPB_SLOTS)
        return;
    elapsed[n] = 0.0;
}

void tacvar_npb_timer_start(int n, int loc_id)
{
    if (n < 0 || n >= TACVAR_NPB_SLOTS)
        return;
    if (!tacvar_is_ready()) {
        if (tacvar_npb_ensure_init() != 0 && !g_init_warned) {
            fprintf(stderr, "tacvar: NPB-OMP init failed; counters/CSV disabled\n");
            g_init_warned = 1;
        }
    }
    note_loc(n, loc_id);
    cpu_start_slot[n] = sched_getcpu();
#if TACVAR_COUNTER_COUNT > 0
    if (tacvar_is_ready())
        TACVAR_COUNTER_READ(counter_start_slot[n]);
#endif
    TACVAR_TIMER_BEGIN(timer_start_raw[n]);
}

void tacvar_npb_timer_stop(int n, int loc_id)
{
    uint64_t timer_stop_raw = 0;
    int64_t elapsed_ns;
    int cpu_stop;
#if TACVAR_COUNTER_COUNT > 0
    uint64_t counter_stop[TACVAR_COUNTER_COUNT];
    uint64_t counter_delta[TACVAR_COUNTER_COUNT];
#endif

    if (n < 0 || n >= TACVAR_NPB_SLOTS)
        return;

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
    note_loc(n, loc_id);

    if (tacvar_is_ready()) {
        tacvar_csv_write_simple(n, loc_id, timer_start_raw[n], timer_stop_raw,
                                elapsed_ns, cpu_start_slot[n], cpu_stop,
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
    if (n < 0 || n >= TACVAR_NPB_SLOTS)
        return 0.0;
    return elapsed[n];
}
