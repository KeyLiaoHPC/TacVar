/**
 * @file counter_papi.c
 * @brief PAPI event-set counter backend for TacVar.
 */
#include "counter_papi.h"
#include "../include/tacvar_measure.h"
#include <papi.h>
#include <stdio.h>
#include <string.h>

#ifndef TACVAR_COUNTER_COUNT
#define TACVAR_COUNTER_COUNT 0
#endif

#if TACVAR_COUNTER_COUNT > 0

static int g_eventset = PAPI_NULL;
static int g_n;
static int g_started;

#ifdef _OPENMP
#include <omp.h>
static int tacvar_papi_thread_id(void)
{
    return omp_get_thread_num();
}
#endif

int tacvar_counter_papi_init(const char *const *names, int n)
{
    int i, code, rc;
    g_n = n;
    g_started = 0;

    rc = PAPI_library_init(PAPI_VER_CURRENT);
    if (rc != PAPI_VER_CURRENT && rc > 0) {
        fprintf(stderr, "tacvar: PAPI_library_init version mismatch\n");
        return -1;
    }
    if (rc < PAPI_OK && rc != PAPI_VER_CURRENT) {
        /* already initialized is OK for some builds */
        if (rc != PAPI_EINVAL) {
            fprintf(stderr, "tacvar: PAPI_library_init failed: %s\n",
                    PAPI_strerror(rc));
            return rc;
        }
    }
#ifdef _OPENMP
    if (PAPI_thread_init((unsigned long (*)(void))tacvar_papi_thread_id) != PAPI_OK) {
        /* non-fatal if already set */
    }
#endif
    g_eventset = PAPI_NULL;
    rc = PAPI_create_eventset(&g_eventset);
    if (rc != PAPI_OK) {
        fprintf(stderr, "tacvar: PAPI_create_eventset: %s\n", PAPI_strerror(rc));
        return rc;
    }
    for (i = 0; i < n; i++) {
        rc = PAPI_event_name_to_code((char *)names[i], &code);
        if (rc != PAPI_OK) {
            fprintf(stderr, "tacvar: PAPI unknown event '%s': %s\n",
                    names[i], PAPI_strerror(rc));
            tacvar_counter_papi_fini();
            return rc;
        }
        rc = PAPI_add_event(g_eventset, code);
        if (rc != PAPI_OK) {
            fprintf(stderr, "tacvar: PAPI_add_event(%s): %s\n",
                    names[i], PAPI_strerror(rc));
            tacvar_counter_papi_fini();
            return rc;
        }
    }
    rc = PAPI_start(g_eventset);
    if (rc != PAPI_OK) {
        fprintf(stderr, "tacvar: PAPI_start: %s\n", PAPI_strerror(rc));
        tacvar_counter_papi_fini();
        return rc;
    }
    g_started = 1;
    return 0;
}

void tacvar_counter_papi_fini(void)
{
    long long discard[TACVAR_MAX_COUNTERS];
    if (g_started && g_eventset != PAPI_NULL) {
        PAPI_stop(g_eventset, discard);
        g_started = 0;
    }
    if (g_eventset != PAPI_NULL) {
        PAPI_cleanup_eventset(g_eventset);
        PAPI_destroy_eventset(&g_eventset);
        g_eventset = PAPI_NULL;
    }
}

void tacvar_counter_papi_read(uint64_t *values)
{
    long long tmp[TACVAR_MAX_COUNTERS];
    int i;
    memset(tmp, 0, sizeof(tmp));
    if (g_eventset != PAPI_NULL)
        PAPI_read(g_eventset, tmp);
    for (i = 0; i < g_n; i++)
        values[i] = (uint64_t)tmp[i];
}

#else

int tacvar_counter_papi_init(const char *const *names, int n)
{
    (void)names; (void)n; return 0;
}
void tacvar_counter_papi_fini(void) {}
void tacvar_counter_papi_read(uint64_t *values) { (void)values; }

#endif
