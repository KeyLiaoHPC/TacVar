/**
 * @file tacvar_csv.c
 * @brief Timestamped data directory, deferred event buffer, and region_info.csv.
 *
 * Event rows are stored in memory during measurement and written only from
 * tacvar_csv_close() / tacvar_fini(). The CSV path is computed at open time
 * but the file itself is created at flush.
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
#include <limits.h>
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

    /* Path only — file is created during final flush. */
    snprintf(st->csv_path, sizeof(st->csv_path),
             "%s/%s_%s_%s_r%04d_t%04d_p%d.csv",
             st->data_dir, suite, bench, klass,
             st->ctx.rank, st->ctx.thread, (int)getpid());
    st->csv_fp = NULL;
    st->header_written = 0;
    st->seq = 0;
    return 0;
}

int tacvar_events_reserve(size_t required_events)
{
    tacvar_event_t *neu;

    if (required_events <= g_tacvar.event_capacity)
        return 0;
    if (required_events > SIZE_MAX / sizeof(tacvar_event_t))
        return -EOVERFLOW;
    neu = (tacvar_event_t *)realloc(g_tacvar.events,
                                    required_events * sizeof(tacvar_event_t));
    if (!neu)
        return -ENOMEM;
    g_tacvar.events = neu;
    g_tacvar.event_capacity = required_events;
    return 0;
}

static int grow_events_unlocked(void)
{
    size_t neu_cap;
    int rc;

    if (g_tacvar.event_capacity == 0)
        neu_cap = 1024;
    else if (g_tacvar.event_capacity > SIZE_MAX / 2)
        return -EOVERFLOW;
    else
        neu_cap = g_tacvar.event_capacity * 2;
    rc = tacvar_events_reserve(neu_cap);
    if (rc != 0 && !g_tacvar.recording_error) {
        fprintf(stderr, "tacvar: event buffer growth failed (%d)\n", rc);
        g_tacvar.recording_error = 1;
    }
    return rc;
}

static void write_header(FILE *fp, tacvar_state_t *st)
{
    int i;
    fprintf(fp,
            "seq,suite,benchmark,class,test_tag,region_id,timer,"
            "raw_start,raw_stop,elapsed_ns,rank,thread,pid,cpu_start,cpu_stop,"
            "migrated,valid,source");
#if TACVAR_COUNTER_COUNT > 0
    fprintf(fp, ",counter_backend");
    for (i = 0; i < st->n_counters; i++) {
        fprintf(fp, ",%s_start,%s_stop,%s_delta",
                st->counter_names[i], st->counter_names[i], st->counter_names[i]);
    }
#else
    (void)i;
    (void)st;
#endif
    fputc('\n', fp);
}

static void write_event_row(FILE *fp, tacvar_state_t *st,
                            const tacvar_event_t *e, uint64_t seq)
{
    int migrated, valid, i;

    migrated = (e->cpu_start >= 0 && e->cpu_stop >= 0 && e->cpu_start != e->cpu_stop);
    valid = !migrated;
    fprintf(fp,
            "%llu,%s,%s,%s,%s,%d,%s,"
            "%llu,%llu,%lld,%d,%d,%d,%d,%d,%d,%d,%s",
            (unsigned long long)seq,
            st->ctx.suite ? st->ctx.suite : "",
            st->ctx.benchmark ? st->ctx.benchmark : "",
            st->ctx.klass ? st->ctx.klass : "",
            st->ctx.test_tag ? st->ctx.test_tag : "",
            e->region_id,
            st->timer_name,
            (unsigned long long)e->raw_start,
            (unsigned long long)e->raw_stop,
            (long long)e->elapsed_ns,
            st->ctx.rank,
            e->thread,
            (int)getpid(),
            e->cpu_start, e->cpu_stop, migrated, valid,
            e->source);
#if TACVAR_COUNTER_COUNT > 0
    fprintf(fp, ",%s", st->counter_backend);
    for (i = 0; i < st->n_counters; i++) {
        fprintf(fp, ",%llu,%llu,%llu",
                (unsigned long long)e->counter_start[i],
                (unsigned long long)e->counter_stop[i],
                (unsigned long long)e->counter_delta[i]);
    }
#else
    (void)i;
#endif
    fputc('\n', fp);
}

static int tacvar_csv_flush(tacvar_state_t *st)
{
    size_t i;
    FILE *fp;

    if (st->event_count == 0)
        return 0;
    if (!st->csv_path[0])
        return -EINVAL;

    fp = fopen(st->csv_path, "w");
    if (!fp) {
        int err = errno ? errno : EIO;
        fprintf(stderr,
                "tacvar: cannot create event CSV '%s': %s "
                "(%zu events still buffered, not discarded)\n",
                st->csv_path, strerror(err), st->event_count);
        return -err;
    }
    write_header(fp, st);
    for (i = 0; i < st->event_count; i++)
        write_event_row(fp, st, &st->events[i], (uint64_t)(i + 1));
    fflush(fp);
    fclose(fp);
    st->csv_fp = NULL;
    st->header_written = 1;
    st->seq = st->event_count;
    return 0;
}

int tacvar_csv_close(tacvar_state_t *st)
{
    int rc;

    if (st->recording_error) {
        fprintf(stderr,
                "tacvar: recording incomplete — events dropped after buffer "
                "growth failure; flushing %zu buffered event(s) "
                "(capacity %zu)\n",
                st->event_count, st->event_capacity);
    }

    rc = tacvar_csv_flush(st);
    if (rc != 0) {
        /* Keep buffer so data is not actively discarded on I/O failure. */
        return rc;
    }

    free(st->events);
    st->events = NULL;
    st->event_count = 0;
    st->event_capacity = 0;
    st->recording_error = 0;
    if (st->csv_fp) {
        fclose(st->csv_fp);
        st->csv_fp = NULL;
    }
    return 0;
}

void tacvar_csv_write_simple(int region_id,
                             uint64_t raw_start, uint64_t raw_stop,
                             int64_t elapsed_ns,
                             int cpu_start, int cpu_stop,
                             const uint64_t *counter_start,
                             const uint64_t *counter_stop,
                             const uint64_t *counter_delta,
                             const char *source)
{
    tacvar_event_t *e;
    int thread_id = 0;
    int i;

    if (!tacvar_is_ready() || g_tacvar.recording_error)
        return;

#ifdef _OPENMP
    thread_id = omp_get_thread_num();
    g_tacvar.ctx.thread = thread_id;
#pragma omp critical(tacvar_csv_write)
#endif
    {
        if (g_tacvar.event_count == g_tacvar.event_capacity) {
            if (grow_events_unlocked() != 0)
                goto out;
        }
        e = &g_tacvar.events[g_tacvar.event_count++];
        memset(e, 0, sizeof(*e));
        e->region_id = region_id;
        e->thread = thread_id;
        e->raw_start = raw_start;
        e->raw_stop = raw_stop;
        e->elapsed_ns = elapsed_ns;
        e->cpu_start = cpu_start;
        e->cpu_stop = cpu_stop;
        if (source)
            snprintf(e->source, sizeof(e->source), "%s", source);
#if TACVAR_COUNTER_COUNT > 0
        for (i = 0; i < g_tacvar.n_counters; i++) {
            e->counter_start[i] = counter_start ? counter_start[i] : 0;
            e->counter_stop[i] = counter_stop ? counter_stop[i] : 0;
            e->counter_delta[i] = counter_delta ? counter_delta[i]
                : (e->counter_stop[i] - e->counter_start[i]);
        }
#else
        (void)counter_start;
        (void)counter_stop;
        (void)counter_delta;
        (void)i;
#endif
out:
        ;
    }
}

void tacvar_region_info_clear(void)
{
    g_tacvar.region_count = 0;
    g_tacvar.region_info_written = 0;
    memset(g_tacvar.regions, 0, sizeof(g_tacvar.regions));
}

void tacvar_region_info_register(const tacvar_region_info_t *regions,
                                 size_t count)
{
    size_t i;
    if (!regions || count == 0)
        return;
    for (i = 0; i < count; i++) {
        if (g_tacvar.region_count >= TACVAR_REGION_INFO_MAX)
            break;
        g_tacvar.regions[g_tacvar.region_count++] = regions[i];
    }
}

static void csv_escape(FILE *fp, const char *s)
{
    const char *p;
    int need_quote = 0;
    if (!s) {
        return;
    }
    for (p = s; *p; p++) {
        if (*p == ',' || *p == '"' || *p == '\n' || *p == '\r') {
            need_quote = 1;
            break;
        }
    }
    if (!need_quote) {
        fputs(s, fp);
        return;
    }
    fputc('"', fp);
    for (p = s; *p; p++) {
        if (*p == '"')
            fputc('"', fp);
        fputc(*p, fp);
    }
    fputc('"', fp);
}

int tacvar_region_info_write(void)
{
    char path[TACVAR_PATH_MAX];
    FILE *fp;
    size_t i;
    struct stat stbuf;

    if (!tacvar_is_ready())
        return -EINVAL;
    if (g_tacvar.ctx.rank != 0 || g_tacvar.ctx.thread != 0)
        return 0;
    if (g_tacvar.region_count == 0)
        return 0;

    snprintf(path, sizeof(path), "%s/region_info.csv", g_tacvar.data_dir);
    if (stat(path, &stbuf) == 0) {
        g_tacvar.region_info_written = 1;
        return 0; /* first writer wins */
    }

    fp = fopen(path, "w");
    if (!fp)
        return -errno;

    fprintf(fp,
            "suite,benchmark,class,test_tag,region_id,region_name,"
            "source_locations,description,active_when,recorded_by_tacvar\n");
    for (i = 0; i < g_tacvar.region_count; i++) {
        const tacvar_region_info_t *r = &g_tacvar.regions[i];
        csv_escape(fp, g_tacvar.ctx.suite);
        fputc(',', fp);
        csv_escape(fp, g_tacvar.ctx.benchmark);
        fputc(',', fp);
        csv_escape(fp, g_tacvar.ctx.klass);
        fputc(',', fp);
        csv_escape(fp, g_tacvar.ctx.test_tag);
        fprintf(fp, ",%d,", r->region_id);
        csv_escape(fp, r->region_name);
        fputc(',', fp);
        csv_escape(fp, r->source_locations);
        fputc(',', fp);
        csv_escape(fp, r->description);
        fputc(',', fp);
        csv_escape(fp, r->active_when);
        fprintf(fp, ",%d\n", r->recorded_by_tacvar);
    }
    fflush(fp);
    fclose(fp);
    g_tacvar.region_info_written = 1;
    return 0;
}
