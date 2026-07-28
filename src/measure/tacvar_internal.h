/**
 * @file tacvar_internal.h
 * @brief Internal state and helpers for TacVar measurement (non-hot-path).
 */
#ifndef TACVAR_INTERNAL_H
#define TACVAR_INTERNAL_H

#include "include/tacvar_measure.h"
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>

#define TACVAR_DIR_NAME_MAX 512
#define TACVAR_PATH_MAX     1024
#define TACVAR_NAME_MAX     128
#define TACVAR_SOURCE_MAX   192
#define TACVAR_REGION_INFO_MAX 64

typedef struct tacvar_event {
    int region_id;
    int thread;
    uint64_t raw_start;
    uint64_t raw_stop;
    int64_t elapsed_ns;
    int cpu_start;
    int cpu_stop;
    uint64_t counter_start[TACVAR_MAX_COUNTERS];
    uint64_t counter_stop[TACVAR_MAX_COUNTERS];
    uint64_t counter_delta[TACVAR_MAX_COUNTERS];
    char source[TACVAR_SOURCE_MAX];
} tacvar_event_t;

typedef struct tacvar_state {
    int initialized;
    pid_t owner_pid;
    tacvar_context_t ctx;
    char data_dir[TACVAR_DIR_NAME_MAX];
    char csv_path[TACVAR_PATH_MAX];
    FILE *csv_fp;                 /**< opened only during final flush */
    int header_written;
    uint64_t seq;
    char timer_name[TACVAR_NAME_MAX];
    char counter_backend[TACVAR_NAME_MAX];
    int n_counters;
    char counter_names[TACVAR_MAX_COUNTERS][TACVAR_NAME_MAX];
    tacvar_event_t *events;
    size_t event_count;
    size_t event_capacity;
    int recording_error;
    tacvar_region_info_t regions[TACVAR_REGION_INFO_MAX];
    size_t region_count;
    int region_info_written;
} tacvar_state_t;

extern tacvar_state_t g_tacvar;

int  tacvar_csv_open(tacvar_state_t *st);
/** Flush then free buffer on success; on flush failure keep buffer and return <0. */
int  tacvar_csv_close(tacvar_state_t *st);
int  tacvar_prepare_data_dir(tacvar_state_t *st, const char *output_root);

uint64_t tacvar_modular_delta(uint64_t start, uint64_t stop, unsigned width_bits);

#endif /* TACVAR_INTERNAL_H */
