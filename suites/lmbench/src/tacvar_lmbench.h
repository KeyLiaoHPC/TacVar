/**
 * @file tacvar_lmbench.h
 * @brief lmbench TacVar sidecar for start/stop.
 */
#ifndef TACVAR_LMBENCH_H
#define TACVAR_LMBENCH_H

#include <sys/time.h>
#include <stdint.h>

void tacvar_lmbench_ensure_init(void);
void tacvar_lmbench_start(struct timeval *tv);
uint64_t tacvar_lmbench_stop(struct timeval *begin, struct timeval *end);
uint64_t tacvar_lmbench_now_us(void);

#endif
