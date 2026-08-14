/**
 * @file tacvar_npb.c
 * @brief NPB-MPI TacVar timers + MPI rank/dir context + timer_info.csv.
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

#define TACVAR_NPB_SLOTS 64

static double elapsed[TACVAR_NPB_SLOTS];
static uint64_t timer_start_raw[TACVAR_NPB_SLOTS];
static int cpu_start_slot[TACVAR_NPB_SLOTS];
static int max_loc[TACVAR_NPB_SLOTS];
#if TACVAR_COUNTER_COUNT > 0
static uint64_t counter_start_slot[TACVAR_NPB_SLOTS][TACVAR_COUNTER_COUNT];
#endif

static char g_bench[64];
static char g_class[8];
static char g_tag[16];
static int g_identity_set;
static int g_init_warned;
static int g_rank;
static int g_timer_info_registered;

typedef struct {
    int region_id;
    const char *name;
} tacvar_npb_name_t;

/* Names for timed slots only (no derived totcomp/totcomm). */
static const tacvar_npb_name_t names_is[] = {
    {0, "total"}, {1, "rcomp"}, {2, "rcomm"}, {3, "verify"},
};
static const tacvar_npb_name_t names_cg[] = {
    {1, "total"}, {2, "conjg"}, {3, "rcomm"}, {4, "ncomm"},
};
static const tacvar_npb_name_t names_mg[] = {
    {1, "total"}, {2, "init"}, {3, "psinv"}, {4, "resid"},
    {5, "rprj3"}, {6, "interp"}, {7, "norm2u3"}, {8, "comm3"}, {9, "rcomm"},
};
static const tacvar_npb_name_t names_ep[] = {
    {1, "total"}, {2, "gpairs"}, {3, "randn"}, {4, "rcomm"},
};
static const tacvar_npb_name_t names_bt[] = {
    {1, "total"}, {2, "i/o"}, {3, "rhs"}, {4, "xsolve"}, {5, "ysolve"},
    {6, "zsolve"}, {7, "bpack"}, {8, "exch"}, {9, "xcomm"}, {10, "ycomm"},
    {11, "zcomm"}, {12, "enorm"}, {13, "iov"},
};
static const tacvar_npb_name_t names_sp[] = {
    {1, "total"}, {2, "rhs"}, {3, "xsolve"}, {4, "ysolve"}, {5, "zsolve"},
    {6, "bpack"}, {7, "exch"}, {8, "xcomm"}, {9, "ycomm"}, {10, "zcomm"},
};
static const tacvar_npb_name_t names_lu[] = {
    {1, "total"}, {2, "rhs"}, {3, "blts"}, {4, "buts"}, {5, "jacld"},
    {6, "jacu"}, {7, "exch"}, {8, "lcomm"}, {9, "ucomm"}, {10, "rcomm"},
};
static const tacvar_npb_name_t names_ft[] = {
    {1, "total"}, {2, "setup"}, {3, "fft"}, {4, "evolve"}, {5, "checksum"},
    {6, "fftlow"}, {7, "fftcopy"}, {8, "transpose"},
    {9, "transxzloc"}, {10, "transxzglo"}, {11, "transxzfin"},
    {12, "transxyloc"}, {13, "transxyglo"}, {14, "transxyfin"},
    {15, "synch"}, {16, "init"},
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
    *n_out = 0;
    return NULL;
}

static void note_loc(int n, int loc_id)
{
    if (n < 0 || n >= TACVAR_NPB_SLOTS)
        return;
    if (loc_id > max_loc[n])
        max_loc[n] = loc_id;
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

    if (g_rank != 0 || !tacvar_is_ready())
        return;

    memset(have_name, 0, sizeof(have_name));
    names = lookup_names(&nnames);
    for (i = 0; i < nnames; i++) {
        int rid = names[i].region_id;
        if (rid >= 0 && rid < TACVAR_NPB_SLOTS) {
            have_name[rid] = 1;
            if (max_loc[rid] < 1)
                max_loc[rid] = 1; /* region exists in catalog */
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
    char path[512], base[256], *b, *p, *rest, *dot2;
    ssize_t n;
    size_t ntag;
    if (g_identity_set)
        return;
    snprintf(g_bench, sizeof(g_bench), "unknown");
    snprintf(g_class, sizeof(g_class), "U");
    g_tag[0] = '\0';
    n = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (n > 0) {
        path[n] = '\0';
        snprintf(base, sizeof(base), "%s", basename(path));
        b = base;
        p = strchr(b, '.');
        if (p) {
            *p = '\0';
            snprintf(g_bench, sizeof(g_bench), "%s", b);
            rest = p + 1;
            if (rest[0]) {
                g_class[0] = rest[0];
                g_class[1] = '\0';
            }
            /* cg.C.x */
            if (rest[0] && rest[1] == '_') {
                dot2 = strchr(rest + 2, '.');
                if (dot2) {
                    ntag = (size_t)(dot2 - (rest + 2));
                    if (ntag >= sizeof(g_tag))
                        ntag = sizeof(g_tag) - 1;
                    memcpy(g_tag, rest + 2, ntag);
                    g_tag[ntag] = '\0';
                }
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
    g_rank = rank;

    memset(&ctx, 0, sizeof(ctx));
    dirbuf[0] = '\0';
    ctx.suite = "npb-mpi";
    ctx.benchmark = g_bench;
    ctx.klass = g_class;
    ctx.test_tag = g_tag[0] ? g_tag : NULL;
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
    if (!g_timer_info_registered) {
        atexit(write_timer_info_atexit);
        g_timer_info_registered = 1;
    }
    return 0;
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
            fprintf(stderr, "tacvar: NPB-MPI init failed; counters/CSV disabled\n");
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
    int64_t elapsed_ns;
    uint64_t timer_stop_raw = 0;
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
