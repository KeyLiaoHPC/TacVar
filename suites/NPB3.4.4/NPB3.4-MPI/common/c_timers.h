/**
 * @file c_timers.h
 * @brief NPB-MPI C timer ABI plus TacVar step/prepare helpers.
 */
#ifndef __C_TIMERS_H
#define __C_TIMERS_H

#include <stdint.h>

extern void   timer_clear( int n );
extern void   timer_start( int n );
extern void   timer_stop( int n );
extern double timer_read( int n );
extern int    check_timer_flag( void );

/* TacVar helpers — instrumentation only; do not alter NPB workload. */
extern void   tacvar_npb_prepare( uint64_t expected_steps );
extern void   tacvar_npb_step_start( void );
extern void   tacvar_npb_step_stop( void );
extern void   tacvar_npb_kernel_start( void );
extern void   tacvar_npb_kernel_stop( void );
extern void   tacvar_npb_set_test_tag( const char *tag );

#endif
