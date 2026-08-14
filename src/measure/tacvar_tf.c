/**
 * @file tacvar_tf.c
 * @brief In-situ TF runtime: ngauge load, per-core ns/tick offset, elapsed_ns.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <sched.h>

#ifdef TACVAR_HAS_GENERATED_CONFIG
#include "tacvar_generated_config.h"
#endif

#include "gauge_tf_insitu.h"

#ifndef TACVAR_TF_SAMPLING
#define TACVAR_TF_SAMPLING 0
#endif
#ifndef TACVAR_TF_TIMER_IS_TICK
#define TACVAR_TF_TIMER_IS_TICK 0
#endif
#ifndef TACVAR_TF_REG_ID
#define TACVAR_TF_REG_ID 0
#endif
#ifndef TACVAR_TF_LOC_ID
#define TACVAR_TF_LOC_ID 0
#endif

#define TACVAR_TF_CPU_MAX 4096

static uint64_t g_tf_ngauge;

#if TACVAR_TF_SAMPLING
static int64_t g_tf_offset[TACVAR_TF_CPU_MAX];
static unsigned char g_tf_valid[TACVAR_TF_CPU_MAX];

static int ci_eq(const char *a, const char *b)
{
    if (!a || !b)
        return 0;
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b))
            return 0;
        a++;
        b++;
    }
    return *a == 0 && *b == 0;
}

static int load_nspg(const char *data_dir, double *nspg_out)
{
    char path[TACVAR_PATH_MAX];
    FILE *fp;
    double v = 0.0;

    snprintf(path, sizeof(path), "%s/nspg.txt", data_dir);
    fp = fopen(path, "r");
    if (!fp)
        return -1;
    if (fscanf(fp, "%lf", &v) != 1 || v <= 0.0) {
        fclose(fp);
        return -1;
    }
    fclose(fp);
    *nspg_out = v;
    return 0;
}

/**
 * Parse median.csv (kernel,class,region_id,loc_id,median[,ngauge])
 * or met_stat.csv (…,region_id,loc_id,…,p50).
 */
static int load_ngauge_from_table(const char *path, const char *bench,
                                  const char *klass, uint64_t *ngauge_out,
                                  double *median_out)
{
    FILE *fp;
    char line[1024], header[1024];
    char *tok, *save;
    int col_k = -1, col_c = -1, col_r = -1, col_l = -1;
    int col_med = -1, col_ng = -1, col_p50 = -1;
    int i, found = 0;
    double med = 0.0;
    uint64_t ng = 0;

    fp = fopen(path, "r");
    if (!fp)
        return -1;
    if (!fgets(header, sizeof(header), fp)) {
        fclose(fp);
        return -1;
    }
    save = NULL;
    i = 0;
    tok = strtok_r(header, ",\r\n", &save);
    while (tok) {
        if (!strcmp(tok, "kernel") || !strcmp(tok, "benchmark"))
            col_k = i;
        else if (!strcmp(tok, "class"))
            col_c = i;
        else if (!strcmp(tok, "region_id"))
            col_r = i;
        else if (!strcmp(tok, "loc_id"))
            col_l = i;
        else if (!strcmp(tok, "median"))
            col_med = i;
        else if (!strcmp(tok, "ngauge"))
            col_ng = i;
        else if (!strcmp(tok, "p50"))
            col_p50 = i;
        i++;
        tok = strtok_r(NULL, ",\r\n", &save);
    }
    if (col_r < 0 || col_l < 0) {
        fclose(fp);
        return -1;
    }
    if (col_med < 0)
        col_med = col_p50;
    if (col_med < 0 && col_ng < 0) {
        fclose(fp);
        return -1;
    }

    while (fgets(line, sizeof(line), fp)) {
        char cols[16][128];
        int n = 0;
        char *p, *sp = NULL;

        if (line[0] == '\0' || line[0] == '\n')
            continue;
        p = strtok_r(line, ",\r\n", &sp);
        while (p && n < 16) {
            snprintf(cols[n], sizeof(cols[n]), "%s", p);
            n++;
            p = strtok_r(NULL, ",\r\n", &sp);
        }
        if (n <= col_r || n <= col_l)
            continue;
        if (atoi(cols[col_r]) != TACVAR_TF_REG_ID)
            continue;
        if (atoi(cols[col_l]) != TACVAR_TF_LOC_ID)
            continue;
        if (col_k >= 0 && n > col_k && bench && bench[0] &&
            !ci_eq(cols[col_k], bench))
            continue;
        if (col_c >= 0 && n > col_c && klass && klass[0] &&
            !ci_eq(cols[col_c], klass))
            continue;
        if (col_med >= 0 && n > col_med)
            med = atof(cols[col_med]);
        if (col_ng >= 0 && n > col_ng)
            ng = (uint64_t)strtoull(cols[col_ng], NULL, 10);
        found = 1;
        break;
    }
    fclose(fp);
    if (!found)
        return -1;
    *median_out = med;
    *ngauge_out = ng;
    return 0;
}

static int load_ngauge(const tacvar_state_t *st)
{
    char path[TACVAR_PATH_MAX];
    const char *bench = st->ctx.benchmark ? st->ctx.benchmark : "";
    const char *klass = st->ctx.klass ? st->ctx.klass : "";
    uint64_t ng = 0;
    double med = 0.0, nspg = 0.0;
    int got;

    snprintf(path, sizeof(path), "%s/median.csv", st->data_dir);
    got = load_ngauge_from_table(path, bench, klass, &ng, &med);
    if (got != 0) {
        snprintf(path, sizeof(path), "%s/met_stat.csv", st->data_dir);
        got = load_ngauge_from_table(path, bench, klass, &ng, &med);
    }
    if (got != 0)
        return -1;
    if (ng >= 1) {
        g_tf_ngauge = ng;
        return 0;
    }
    if (load_nspg(st->data_dir, &nspg) != 0 || nspg <= 0.0 || med <= 0.0)
        return -1;
    ng = (uint64_t)(med / nspg + 0.5);
    if (ng < 1)
        return -1;
    g_tf_ngauge = ng;
    return 0;
}
#endif /* TACVAR_TF_SAMPLING */

int tacvar_tf_prepare(void)
{
#if !TACVAR_TF_SAMPLING
    g_tf_ngauge = 0;
    return 0;
#else
    int rc;

    memset(g_tf_valid, 0, sizeof(g_tf_valid));
    g_tf_ngauge = 0;

#ifdef __x86_64__
    rc = tacvar_tsc_calibrate();
    if (rc != 0) {
        fprintf(stderr, "tacvar: TSC calibrate failed for TF tick-to-ns (%d)\n", rc);
        return rc;
    }
#endif
    rc = load_ngauge(&g_tacvar);
    if (rc != 0) {
        fprintf(stderr,
                "tacvar: TF ngauge missing (need %s/median.csv "
                "region_id=%d loc_id=%d and nspg.txt)\n",
                g_tacvar.data_dir, TACVAR_TF_REG_ID, TACVAR_TF_LOC_ID);
        return -1;
    }
    return 0;
#endif
}

uint64_t tacvar_tf_ngauge(void)
{
    return g_tf_ngauge;
}

void tacvar_tf_ensure_offset(int cpu)
{
#if !TACVAR_TF_SAMPLING || TACVAR_TF_TIMER_IS_TICK
    (void)cpu;
    return;
#else
    uint64_t o = 0, k = 0;

    if (cpu < 0 || cpu >= TACVAR_TF_CPU_MAX)
        return;
    if (g_tf_valid[cpu])
        return;
    TACVAR_TIMER_BEGIN(o);
    TACVAR_TF_TICK_READ(k);
    g_tf_offset[cpu] = (int64_t)o - TACVAR_TF_TICK_TO_NS(k);
    g_tf_valid[cpu] = 1;
#endif
}

int64_t tacvar_tf_elapsed_ns(int orig_is_start, uint64_t t_orig, uint64_t t_tick)
{
#if !TACVAR_TF_SAMPLING
    (void)orig_is_start;
    (void)t_orig;
    (void)t_tick;
    return 0;
#elif TACVAR_TF_TIMER_IS_TICK
    uint64_t dt;

    (void)orig_is_start;
    dt = (t_orig > t_tick) ? (t_orig - t_tick) : (t_tick - t_orig);
    return TACVAR_TF_TICK_TO_NS(dt);
#else
    int cpu = sched_getcpu();
    int64_t tick_ns, orig_ns, el;

    if (cpu < 0 || cpu >= TACVAR_TF_CPU_MAX || !g_tf_valid[cpu])
        return 0;
    tick_ns = TACVAR_TF_TICK_TO_NS(t_tick) + g_tf_offset[cpu];
    orig_ns = (int64_t)t_orig;
    el = orig_is_start ? (tick_ns - orig_ns) : (orig_ns - tick_ns);
    return el < 0 ? 0 : el;
#endif
}
