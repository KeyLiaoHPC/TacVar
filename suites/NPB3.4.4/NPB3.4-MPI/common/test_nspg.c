/**
 * @file test_nspg.c
 * @brief Fit nanoseconds per subtraction-gauge step using TACVAR_TIMER macros.
 *
 * NPB-MPI copy under common/; build bin/test_nspg.x with `make nspg`.
 * Doubling nsub from 1000; 1s warmup then min of NTEST samples per step.
 * From the 5th group, OLS slope of the latest 5 (nsub, tmin) points is nspg.
 * Stop when R^2 >= R2_THRS (default 0.999999). Rank 0 writes nspg.txt.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
#include "tacvar_generated_config.h"
#include "gauge_sub.h"

#ifndef NTEST
#define NTEST 10000
#endif
#ifndef TMAX_NS
#define TMAX_NS 100000000LL
#endif
#ifndef NSTEPS_MAX
#define NSTEPS_MAX 32
#endif
#ifndef NFIT
#define NFIT 5
#endif
#ifndef NG_START
#define NG_START 1000ULL
#endif
#ifndef WARMUP_NS
#define WARMUP_NS 1000000000LL
#endif
#ifndef R2_THRS
#define R2_THRS 0.999999
#endif

static int64_t run_once(uint64_t nsub)
{
    uint64_t t0 = 0, t1 = 0;
    TACVAR_TIMER_BEGIN(t0);
    TACVAR_TF_GAUGE_N(nsub);
    TACVAR_TIMER_END(t1);
    return TACVAR_TIMER_DELTA_NS(t0, t1);
}

static void warmup_1s(uint64_t nsub)
{
    struct timespec t0, t1;
    int64_t elapsed;

    clock_gettime(CLOCK_MONOTONIC, &t0);
    do {
        (void)run_once(nsub);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        elapsed = (int64_t)(t1.tv_sec - t0.tv_sec) * 1000000000LL
                + (int64_t)(t1.tv_nsec - t0.tv_nsec);
    } while (elapsed < WARMUP_NS);
}

/**
 * OLS y = a*x + b on the latest NFIT points. nspg = slope a.
 * Returns 0 on success, -1 if the x-denominator is 0.
 */
static int fit_ols_window(const uint64_t *p_ng, const int64_t *p_tm_min,
                          int n, double *nspg_out, double *r2_out)
{
    int i, i0 = n - NFIT;
    double sx = 0.0, sy = 0.0, sxy = 0.0, sx2 = 0.0;
    double k = (double)NFIT, denom, a, b, mean_y, ss_res = 0.0, ss_tot = 0.0;

    for (i = i0; i < n; i++) {
        double x = (double)p_ng[i];
        double y = (double)p_tm_min[i];
        sx += x;
        sy += y;
        sxy += x * y;
        sx2 += x * x;
    }
    denom = k * sx2 - sx * sx;
    if (denom == 0.0)
        return -1;
    a = (k * sxy - sx * sy) / denom;
    b = (sy - a * sx) / k;
    mean_y = sy / k;
    for (i = i0; i < n; i++) {
        double y = (double)p_tm_min[i];
        double yhat = a * (double)p_ng[i] + b;
        ss_res += (y - yhat) * (y - yhat);
        ss_tot += (y - mean_y) * (y - mean_y);
    }
    *nspg_out = a;
    if (ss_tot == 0.0)
        *r2_out = 1.0;
    else
        *r2_out = 1.0 - ss_res / ss_tot;
    return 0;
}

int main(int argc, char **argv)
{
    int rank = 0, nrank = 1, rc, have_nspg = 0, fit_ok, r2_met = 0;
    uint64_t ng;
    int n = 0, i;
    int64_t tpre, tmin, sample;
    int64_t p_tm_min[NSTEPS_MAX];
    uint64_t p_ng[NSTEPS_MAX];
    double nspg = 0.0, rsquare = 0.0;
    const char *out_path = "nspg.txt";
    FILE *fp;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    if (argc > 1 && argv[1][0])
        out_path = argv[1];

    rc = TACVAR_TIMER_INIT();
    if (rc != 0) {
        if (rank == 0)
            fprintf(stderr, "test_nspg: timer init failed (%d)\n", rc);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    ng = NG_START;
    tpre = 0;
    n = 0;
    while (tpre <= TMAX_NS && n < NSTEPS_MAX) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        warmup_1s(ng);
        tmin = 0;
        for (i = 0; i < NTEST; i++) {
            sample = run_once(ng);
            if (i == 0 || sample < tmin)
                tmin = sample;
        }
        MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_INT64_T, MPI_MIN,
                      MPI_COMM_WORLD);
        p_tm_min[n] = tmin;
        p_ng[n] = ng;
        n++;
        fit_ok = 0;
        if (n >= NFIT) {
            if (fit_ols_window(p_ng, p_tm_min, n, &nspg, &rsquare) == 0) {
                have_nspg = 1;
                fit_ok = 1;
            } else if (rank == 0) {
                fprintf(stderr,
                        "test_nspg: OLS denominator is 0, skip nspg update"
                        " (nsub=%" PRIu64 ")\n", ng);
            }
        }
        if (rank == 0) {
            if (fit_ok)
                printf("nsub=%" PRIu64 " tmin=%" PRIi64
                       " ns  nfit=%d  nspg=%.9f  R^2=%.9f\n",
                       ng, tmin, NFIT, nspg, rsquare);
            else
                printf("nsub=%" PRIu64 " tmin=%" PRIi64 " ns\n", ng, tmin);
            fflush(stdout);
        }
        if (fit_ok && rsquare >= R2_THRS) {
            r2_met = 1;
            break;
        }
        ng *= 2;
        tpre = tmin * 2;
        if (tmin <= 0)
            break;
    }

    if (n < NFIT || !have_nspg || nspg <= 0.0) {
        if (rank == 0) {
            if (n < NFIT)
                fprintf(stderr,
                        "test_nspg: need at least %d groups to fit nspg (got %d)\n",
                        NFIT, n);
            else
                fprintf(stderr, "test_nspg: failed to fit nspg\n");
        }
        TACVAR_TIMER_FINI();
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        if (!r2_met)
            fprintf(stderr,
                    "test_nspg: R^2 did not reach %.9f (last R^2=%.9f); "
                    "using last nspg\n",
                    R2_THRS, rsquare);
        printf("nspg=%.9f ns/step (nrank=%d nsteps=%d R^2=%.9f)\n",
               nspg, nrank, n, rsquare);
        fp = fopen(out_path, "w");
        if (!fp) {
            perror(out_path);
            TACVAR_TIMER_FINI();
            MPI_Finalize();
            return 1;
        }
        fprintf(fp, "%.12g\n", nspg);
        fclose(fp);
        printf("wrote %s\n", out_path);
    }

    TACVAR_TIMER_FINI();
    MPI_Finalize();
    return 0;
}
