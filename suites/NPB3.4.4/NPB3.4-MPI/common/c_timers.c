/**
 * @file c_timers.c
 * @brief NPB-MPI C timer ABI — delegates to tacvar_npb_*.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tacvar_npb.h"

void timer_clear(int n) { tacvar_npb_timer_clear(n); }
void timer_start(int n, int loc_id) { tacvar_npb_timer_start(n, loc_id); }
void timer_stop(int n, int loc_id)  { tacvar_npb_timer_stop(n, loc_id); }
double timer_read(int n) { return tacvar_npb_timer_read(n); }

int check_timer_flag(void)
{
    int timer_on = 0;
    char *ev = getenv("NPB_TIMER_FLAG");

    if (ev) {
        if (*ev == '\0')
            timer_on = 1;
        else if (*ev >= '1' && *ev <= '9')
            timer_on = 1;
        else if (strcmp(ev, "on") == 0 || strcmp(ev, "ON") == 0 ||
                 strcmp(ev, "yes") == 0 || strcmp(ev, "YES") == 0 ||
                 strcmp(ev, "true") == 0 || strcmp(ev, "TRUE") == 0)
            timer_on = 1;
    } else {
        FILE *fp = fopen("timer.flag", "r");
        if (fp != NULL) {
            fclose(fp);
            timer_on = 1;
        }
    }
    return timer_on;
}
