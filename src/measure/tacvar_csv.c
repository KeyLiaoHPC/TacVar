/**
 * @file tacvar_csv.c
 * @brief Timestamped data directory and per-writer CSV output.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "tacvar_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <sched.h>
#ifdef _OPENMP
#include <omp.h>
#endif

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

int tacvar_prepare_data_dir(tacvar_state_t *st, const char *output_root)
{
    const char *env_dir;
    const char *root;
    time_t now;
    struct tm tm_now;
    char stamp[32];
    int rc;

    env_dir = getenv("TACVAR_DATA_DIR");
    if (env_dir && env_dir[0]) {
        snprintf(st->data_dir, sizeof(st->data_dir), "%s", env_dir);
        rc = mkdir(st->data_dir, 0755);
        if (rc != 0 && errno != EEXIST)
            return -errno;
        return 0;
    }

    root = output_root && output_root[0] ? output_root : TACVAR_OUTPUT_ROOT_DEFAULT;
    now = time(NULL);
    localtime_r(&now, &tm_now);
    strftime(stamp, sizeof(stamp), "%Y%m%dT%H%M%S", &tm_now);
    snprintf(st->data_dir, sizeof(st->data_dir), "%s/data_%s", root, stamp);
    rc = mkdir(st->data_dir, 0755);
    if (rc != 0 && errno != EEXIST) {
        /* collision: append pid */
        snprintf(st->data_dir, sizeof(st->data_dir), "%s/data_%s_%d",
                 root, stamp, (int)getpid());
        rc = mkdir(st->data_dir, 0755);
        if (rc != 0 && errno != EEXIST)
            return -errno;
    }
    /* Publish for children after fork/exec. */
    setenv("TACVAR_DATA_DIR", st->data_dir, 1);
    return 0;
}

int tacvar_csv_open(tacvar_state_t *st)
{
    const char *suite = st->ctx.suite ? st->ctx.suite : "suite";
    const char *bench = st->ctx.benchmark ? st->ctx.benchmark : "bench";
    const char *klass = st->ctx.klass ? st->ctx.klass : "X";

    snprintf(st->csv_path, sizeof(st->csv_path),
             "%s/%s_%s_%s_r%04d_t%04d_p%d.csv",
             st->data_dir, suite, bench, klass,
             st->ctx.rank, st->ctx.thread, (int)getpid());
    st->csv_fp = fopen(st->csv_path, "a");
    if (!st->csv_fp)
        return -errno;
    setvbuf(st->csv_fp, NULL, _IOFBF, 1 << 16);
    st->header_written = 0;
    st->seq = 0;
    return 0;
}

void tacvar_csv_close(tacvar_state_t *st)
{
    if (st->csv_fp) {
        fflush(st->csv_fp);
        fclose(st->csv_fp);
        st->csv_fp = NULL;
    }
}

static void write_header(tacvar_state_t *st)
{
    int i;
    fprintf(st->csv_fp,
            "seq,suite,benchmark,class,test_tag,region_id,timer,"
            "raw_start,raw_stop,elapsed_ns,rank,thread,pid,cpu_start,cpu_stop,"
            "migrated,valid");
#if TACVAR_COUNTER_COUNT > 0
    fprintf(st->csv_fp, ",counter_backend");
    for (i = 0; i < st->n_counters; i++) {
        fprintf(st->csv_fp, ",%s_start,%s_stop,%s_delta",
                st->counter_names[i], st->counter_names[i], st->counter_names[i]);
    }
#else
    (void)i;
#endif
    fputc('\n', st->csv_fp);
    st->header_written = 1;
}

void tacvar_csv_write_simple(int region_id,
                             uint64_t raw_start, uint64_t raw_stop,
                             int64_t elapsed_ns,
                             int cpu_start, int cpu_stop,
                             const uint64_t *counter_start,
                             const uint64_t *counter_stop,
                             const uint64_t *counter_delta)
{
    int migrated, valid, i;
    if (!g_tacvar.csv_fp)
        return;

#ifdef _OPENMP
    g_tacvar.ctx.thread = omp_get_thread_num();
#pragma omp critical(tacvar_csv_write)
#endif
    {
    if (!g_tacvar.header_written)
        write_header(&g_tacvar);

    migrated = (cpu_start >= 0 && cpu_stop >= 0 && cpu_start != cpu_stop);
    valid = !migrated;
    g_tacvar.seq++;

    fprintf(g_tacvar.csv_fp,
            "%llu,%s,%s,%s,%s,%d,%s,"
            "%llu,%llu,%lld,%d,%d,%d,%d,%d,%d,%d",
            (unsigned long long)g_tacvar.seq,
            g_tacvar.ctx.suite ? g_tacvar.ctx.suite : "",
            g_tacvar.ctx.benchmark ? g_tacvar.ctx.benchmark : "",
            g_tacvar.ctx.klass ? g_tacvar.ctx.klass : "",
            g_tacvar.ctx.test_tag ? g_tacvar.ctx.test_tag : "",
            region_id,
            g_tacvar.timer_name,
            (unsigned long long)raw_start,
            (unsigned long long)raw_stop,
            (long long)elapsed_ns,
            g_tacvar.ctx.rank,
            g_tacvar.ctx.thread,
            (int)getpid(),
            cpu_start, cpu_stop, migrated, valid);
#if TACVAR_COUNTER_COUNT > 0
    fprintf(g_tacvar.csv_fp, ",%s", g_tacvar.counter_backend);
    for (i = 0; i < g_tacvar.n_counters; i++) {
        uint64_t s = counter_start ? counter_start[i] : 0;
        uint64_t e = counter_stop ? counter_stop[i] : 0;
        uint64_t d = counter_delta ? counter_delta[i] : (e - s);
        fprintf(g_tacvar.csv_fp, ",%llu,%llu,%llu",
                (unsigned long long)s, (unsigned long long)e,
                (unsigned long long)d);
    }
#else
    (void)counter_start;
    (void)counter_stop;
    (void)counter_delta;
    (void)i;
#endif
    fputc('\n', g_tacvar.csv_fp);
    }
}
