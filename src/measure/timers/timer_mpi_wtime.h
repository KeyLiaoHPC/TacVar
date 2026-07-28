/**
 * @file timer_mpi_wtime.h
 * @brief MPI_Wtime timer; raw is nanoseconds. NPB-MPI only.
 */
#ifndef TACVAR_TIMER_MPI_WTIME_H
#define TACVAR_TIMER_MPI_WTIME_H

#include <mpi.h>
#include <stdint.h>

static inline int tacvar_timer_mpi_wtime_init(void) { return 0; }
static inline void tacvar_timer_mpi_wtime_fini(void) {}

static inline void tacvar_timer_mpi_wtime_begin(uint64_t *raw)
{
    *raw = (uint64_t)(MPI_Wtime() * 1e9);
}

static inline void tacvar_timer_mpi_wtime_end(uint64_t *raw)
{
    *raw = (uint64_t)(MPI_Wtime() * 1e9);
}

static inline int64_t tacvar_timer_mpi_wtime_delta_ns(uint64_t b, uint64_t e)
{
    return (int64_t)(e - b);
}

#endif
