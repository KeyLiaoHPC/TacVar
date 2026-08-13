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

typedef struct tacvar_state {
    int initialized;
    pid_t owner_pid;
    tacvar_context_t ctx;
    char data_dir[TACVAR_DIR_NAME_MAX];
    char kernel_dir[TACVAR_DIR_NAME_MAX]; /**< data_dir/Kernel.CLASS */
    char short_host[TACVAR_NAME_MAX];
    char csv_path[TACVAR_PATH_MAX];
    FILE *csv_fp;
    int header_written;
    uint64_t seq;
    char timer_name[TACVAR_NAME_MAX];
    char counter_backend[TACVAR_NAME_MAX];
    int n_counters;
    char counter_names[TACVAR_MAX_COUNTERS][TACVAR_NAME_MAX];
} tacvar_state_t;

extern tacvar_state_t g_tacvar;

int  tacvar_csv_open(tacvar_state_t *st);
void tacvar_csv_close(tacvar_state_t *st);
int  tacvar_prepare_data_dir(tacvar_state_t *st, const char *output_root);
int  tacvar_prepare_kernel_dir(tacvar_state_t *st);

/* Optional weak hooks filled by selected backend .c files. */
int  tacvar_timer_backend_init(void);
void tacvar_timer_backend_fini(void);
int  tacvar_counter_backend_init(const char *const *names, int n);
void tacvar_counter_backend_fini(void);

uint64_t tacvar_modular_delta(uint64_t start, uint64_t stop, unsigned width_bits);

#endif /* TACVAR_INTERNAL_H */
