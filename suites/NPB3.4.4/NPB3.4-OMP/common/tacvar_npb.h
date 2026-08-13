/**
 * @file tacvar_npb.h
 * @brief NPB TacVar timer implementation shared by C and Fortran.
 */
#ifndef TACVAR_NPB_H
#define TACVAR_NPB_H

#ifdef __cplusplus
extern "C" {
#endif

void tacvar_npb_timer_clear(int n);
void tacvar_npb_timer_start(int n, int loc_id);
void tacvar_npb_timer_stop(int n, int loc_id);
double tacvar_npb_timer_read(int n);
int tacvar_npb_ensure_init(void);

#ifdef __cplusplus
}
#endif
#endif
