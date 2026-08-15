/*
 * timer_rate.c — runtime rate (ns per hardware tick) loader.
 *
 * Priority:  caller-provided file > TACVAR_NSPT_FILE macro > fallback.
 * The file format is what probe_counter_freq --nspt-out writes:
 *   # comment lines...
 *   10.000754003
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "timer_rate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double g_tacvar_ns_per_tick = 0.0;

#ifndef TACVAR_NSPT_FILE
#define TACVAR_NSPT_FILE ""
#endif

int tacvar_rate_load_file(const char *nspt_file)
{
    FILE *f;
    char line[256];
    double v = 0.0;

    if (!nspt_file || !*nspt_file)
        return -1;
    f = fopen(nspt_file, "r");
    if (!f)
        return -1;
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        char *end;
        double x;
        /* skip blank and comment lines */
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0')
            continue;
        x = strtod(p, &end);
        if (end != p && x > 0.0) {
            v = x;
            break;
        }
    }
    fclose(f);
    if (v <= 0.0)
        return -1;
    g_tacvar_ns_per_tick = v;
    return 0;
}

int tacvar_rate_init(const char *nspt_file, double fallback)
{
    const char *path = nspt_file;

    if (!path || !*path)
        path = TACVAR_NSPT_FILE;
    if (path && *path) {
        if (tacvar_rate_load_file(path) == 0) {
            fprintf(stderr, "tacvar: ns_per_tick=%.9f (file %s)\n",
                    g_tacvar_ns_per_tick, path);
            return 0;
        }
        fprintf(stderr,
                "tacvar: warn: cannot read nspt file '%s'; "
                "using fallback %.9f ns/tick\n", path, fallback);
    }
    g_tacvar_ns_per_tick = fallback;
    return 0;
}