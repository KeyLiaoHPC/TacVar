/**
 * @file tacvar_csv.c
 * @brief Timestamped data directory and per-writer CSV output.
 *
 * Layout: DATA_ROOT/timer_info.csv
 *         DATA_ROOT/Kernel.CLASS/<short_host>_rRRRR_tTTTT_pPID.csv
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
#include <sys/file.h>
#include <fcntl.h>
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

#ifdef _OPENMP
static FILE *tl_csv_fp;
static char tl_csv_path[TACVAR_PATH_MAX];
static int tl_header_written;
static uint64_t tl_seq;
#pragma omp threadprivate(tl_csv_fp, tl_csv_path, tl_header_written, tl_seq)
#endif

static void fill_short_host(char *buf, size_t buflen)
{
    char host[256];
    char *dot;
    size_t i;

    if (gethostname(host, sizeof(host)) != 0) {
        snprintf(buf, buflen, "unknown");
        return;
    }
    host[sizeof(host) - 1] = '\0';
    dot = strchr(host, '.');
    if (dot)
        *dot = '\0';
    for (i = 0; host[i]; i++) {
        if (host[i] == '/')
            host[i] = '_';
    }
    if (!host[0])
        snprintf(buf, buflen, "unknown");
    else
        snprintf(buf, buflen, "%s", host);
}

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

#if defined(TACVAR_TF_SAMPLING) && TACVAR_TF_SAMPLING
    if (TACVAR_TF_DATA_ROOT[0]) {
        snprintf(st->data_dir, sizeof(st->data_dir), "%s", TACVAR_TF_DATA_ROOT);
        rc = mkdir(st->data_dir, 0755);
        if (rc != 0 && errno != EEXIST)
            return -errno;
        setenv("TACVAR_DATA_DIR", st->data_dir, 1);
        return 0;
    }
#endif

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

int tacvar_prepare_kernel_dir(tacvar_state_t *st)
{
    const char *bench = st->ctx.benchmark ? st->ctx.benchmark : "bench";
    const char *klass = st->ctx.klass ? st->ctx.klass : "X";
    const char *tag = st->ctx.test_tag;
    int rc;

    fill_short_host(st->short_host, sizeof(st->short_host));
    if (tag && tag[0])
        snprintf(st->kernel_dir, sizeof(st->kernel_dir), "%s/%s.%s_%s",
                 st->data_dir, bench, klass, tag);
    else
        snprintf(st->kernel_dir, sizeof(st->kernel_dir), "%s/%s.%s",
                 st->data_dir, bench, klass);
    rc = mkdir(st->kernel_dir, 0755);
    if (rc != 0 && errno != EEXIST)
        return -errno;
    return 0;
}

static int open_csv_at(tacvar_state_t *st, FILE **fp, char *path, size_t pathlen,
                       int *header_written, uint64_t *seq)
{
    int thread = st->ctx.thread;
    int rc;

    rc = tacvar_prepare_kernel_dir(st);
    if (rc != 0)
        return rc;

#ifdef _OPENMP
    thread = omp_get_thread_num();
    st->ctx.thread = thread;
#endif

    snprintf(path, pathlen, "%s/%s_r%04d_t%04d_p%d.csv",
             st->kernel_dir, st->short_host,
             st->ctx.rank, thread, (int)getpid());
    *fp = fopen(path, "a");
    if (!*fp)
        return -errno;
    setvbuf(*fp, NULL, _IOFBF, 1 << 16);
    *header_written = 0;
    *seq = 0;
    return 0;
}

int tacvar_csv_open(tacvar_state_t *st)
{
#ifdef _OPENMP
    /* Defer open until first write so tNNNN matches the real OMP thread. */
    (void)st;
    return tacvar_prepare_kernel_dir(st);
#else
    return open_csv_at(st, &st->csv_fp, st->csv_path, sizeof(st->csv_path),
                       &st->header_written, &st->seq);
#endif
}

void tacvar_csv_close(tacvar_state_t *st)
{
#ifdef _OPENMP
    if (tl_csv_fp) {
        fflush(tl_csv_fp);
        fclose(tl_csv_fp);
        tl_csv_fp = NULL;
    }
#endif
    if (st->csv_fp) {
        fflush(st->csv_fp);
        fclose(st->csv_fp);
        st->csv_fp = NULL;
    }
}

static void write_header(FILE *fp, tacvar_state_t *st)
{
    int i;
#if TACVAR_COUNTER_COUNT <= 0
    (void)st;
#endif
    fprintf(fp,
            "seq,region_id,loc_id,raw_start,raw_stop,elapsed_ns,"
            "rank,thread,pid,cpu_start,cpu_stop,migrated,valid");
#if TACVAR_COUNTER_COUNT > 0
    fprintf(fp, ",counter_backend");
    for (i = 0; i < st->n_counters; i++) {
        fprintf(fp, ",%s_start,%s_stop,%s_delta",
                st->counter_names[i], st->counter_names[i], st->counter_names[i]);
    }
#else
    (void)i;
#endif
    fputc('\n', fp);
}

void tacvar_csv_write_simple(int region_id, int loc_id,
                             uint64_t raw_start, uint64_t raw_stop,
                             int64_t elapsed_ns,
                             int cpu_start, int cpu_stop,
                             const uint64_t *counter_start,
                             const uint64_t *counter_stop,
                             const uint64_t *counter_delta)
{
    int migrated, valid, i;
    FILE *fp;
    int *hdr;
    uint64_t *seq;

    if (!g_tacvar.initialized)
        return;

#ifdef _OPENMP
    g_tacvar.ctx.thread = omp_get_thread_num();
    if (!tl_csv_fp) {
#pragma omp critical(tacvar_csv_mkdir)
        {
            (void)tacvar_prepare_kernel_dir(&g_tacvar);
        }
        if (open_csv_at(&g_tacvar, &tl_csv_fp, tl_csv_path, sizeof(tl_csv_path),
                        &tl_header_written, &tl_seq) != 0)
            return;
    }
    fp = tl_csv_fp;
    hdr = &tl_header_written;
    seq = &tl_seq;
#else
    fp = g_tacvar.csv_fp;
    hdr = &g_tacvar.header_written;
    seq = &g_tacvar.seq;
    if (!fp)
        return;
#endif

    if (!*hdr) {
        write_header(fp, &g_tacvar);
        *hdr = 1;
    }

    migrated = (cpu_start >= 0 && cpu_stop >= 0 && cpu_start != cpu_stop);
    valid = !migrated;
    (*seq)++;

    fprintf(fp,
            "%llu,%d,%d,%llu,%llu,%lld,%d,%d,%d,%d,%d,%d,%d",
            (unsigned long long)*seq,
            region_id,
            loc_id,
            (unsigned long long)raw_start,
            (unsigned long long)raw_stop,
            (long long)elapsed_ns,
            g_tacvar.ctx.rank,
            g_tacvar.ctx.thread,
            (int)getpid(),
            cpu_start, cpu_stop, migrated, valid);
#if TACVAR_COUNTER_COUNT > 0
    fprintf(fp, ",%s", g_tacvar.counter_backend);
    for (i = 0; i < g_tacvar.n_counters; i++) {
        uint64_t s = counter_start ? counter_start[i] : 0;
        uint64_t e = counter_stop ? counter_stop[i] : 0;
        uint64_t d = counter_delta ? counter_delta[i] : (e - s);
        fprintf(fp, ",%llu,%llu,%llu",
                (unsigned long long)s, (unsigned long long)e,
                (unsigned long long)d);
    }
#else
    (void)counter_start;
    (void)counter_stop;
    (void)counter_delta;
    (void)i;
#endif
    fputc('\n', fp);
}

int tacvar_write_timer_info(const int *region_ids,
                            const int *nlocs,
                            const char *const *names,
                            int n)
{
    char path[TACVAR_PATH_MAX];
    FILE *fp;
    int i, fd, rc;
    off_t sz;
    const char *suite;
    const char *bench;
    const char *klass;
    const char *tag;
    const char *timer;

    if (!g_tacvar.initialized || n <= 0 || !region_ids || !nlocs || !names)
        return -EINVAL;
    if (!g_tacvar.data_dir[0])
        return -EINVAL;

    suite = g_tacvar.ctx.suite ? g_tacvar.ctx.suite : "";
    bench = g_tacvar.ctx.benchmark ? g_tacvar.ctx.benchmark : "";
    klass = g_tacvar.ctx.klass ? g_tacvar.ctx.klass : "";
    tag = g_tacvar.ctx.test_tag ? g_tacvar.ctx.test_tag : "";
    timer = g_tacvar.timer_name[0] ? g_tacvar.timer_name : "";

    snprintf(path, sizeof(path), "%s/timer_info.csv", g_tacvar.data_dir);
    fd = open(path, O_RDWR | O_CREAT, 0644);
    if (fd < 0)
        return -errno;
    if (flock(fd, LOCK_EX) != 0) {
        rc = -errno;
        close(fd);
        return rc;
    }
    sz = lseek(fd, 0, SEEK_END);
    if (sz < 0) {
        rc = -errno;
        flock(fd, LOCK_UN);
        close(fd);
        return rc;
    }
    fp = fdopen(fd, "a");
    if (!fp) {
        rc = -errno;
        flock(fd, LOCK_UN);
        close(fd);
        return rc;
    }
    if (sz == 0)
        fprintf(fp, "suite,benchmark,class,test_tag,timer,region_id,nloc,name\n");
    for (i = 0; i < n; i++) {
        fprintf(fp, "%s,%s,%s,%s,%s,%d,%d,%s\n",
                suite, bench, klass, tag, timer,
                region_ids[i], nlocs[i],
                names[i] ? names[i] : "");
    }
    fflush(fp);
    flock(fd, LOCK_UN);
    fclose(fp);
    return 0;
}
