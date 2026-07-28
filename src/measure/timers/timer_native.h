/**
 * @file timer_native.h
 * @brief Suite-native timer alias resolved at generate time via TACVAR_CONSUMER.
 *
 * Includes the original suite timer implementation (compile-time only).
 */
#ifndef TACVAR_TIMER_NATIVE_H
#define TACVAR_TIMER_NATIVE_H

#if !defined(TACVAR_CONSUMER)
#error "TACVAR_CONSUMER required for native timer"
#endif

#if TACVAR_CONSUMER == TACVAR_CONSUMER_NPB_MPI
#include "timer_mpi_wtime.h"
#define tacvar_timer_native_init      tacvar_timer_mpi_wtime_init
#define tacvar_timer_native_fini      tacvar_timer_mpi_wtime_fini
#define tacvar_timer_native_begin     tacvar_timer_mpi_wtime_begin
#define tacvar_timer_native_end       tacvar_timer_mpi_wtime_end
#define tacvar_timer_native_delta_ns  tacvar_timer_mpi_wtime_delta_ns
#elif TACVAR_CONSUMER == TACVAR_CONSUMER_NPB_OMP
#include "timer_omp_get_wtime.h"
#define tacvar_timer_native_init      tacvar_timer_omp_get_wtime_init
#define tacvar_timer_native_fini      tacvar_timer_omp_get_wtime_fini
#define tacvar_timer_native_begin     tacvar_timer_omp_get_wtime_begin
#define tacvar_timer_native_end       tacvar_timer_omp_get_wtime_end
#define tacvar_timer_native_delta_ns  tacvar_timer_omp_get_wtime_delta_ns
#elif TACVAR_CONSUMER == TACVAR_CONSUMER_LMBENCH
#include "timer_gettimeofday.h"
#define tacvar_timer_native_init      tacvar_timer_gettimeofday_init
#define tacvar_timer_native_fini      tacvar_timer_gettimeofday_fini
#define tacvar_timer_native_begin     tacvar_timer_gettimeofday_begin
#define tacvar_timer_native_end       tacvar_timer_gettimeofday_end
#define tacvar_timer_native_delta_ns  tacvar_timer_gettimeofday_delta_ns
#else
#error "Unknown TACVAR_CONSUMER for native timer"
#endif

#endif
