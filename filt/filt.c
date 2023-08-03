/**
 * @file filtv3.c
 * @author Key Liao
 *
 */
#define _GNU_SOURCE
#define _ISOC11_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <error.h>
#include <errno.h>
#include <math.h>
#include <argp.h>
#include <sys/time.h>
#include "filt.h"


/*=== ARGP ===*/
static struct argp_option options[] = {
    {"met-file", 'm', "FILE", 0, "Input measurement file", 0},
    {"sample-file", 's', "FILE", 0, "Input timing fluctuation file", 0},
    {"base-time", 'b', "TIME", 0, "Base time of timing fluctuation (ns)", 0},
    {"step-time", 't', "TIME", 0, "Time step in initial guess (ns)", 0},
    {"round", 'r', "ROUND", 0, "Time for rounding data bin edge (ns)", 0},
    {"nsamp", 'n', "NSAMP", 0, "Number of samples in each optimization step", 0},
    {"plow", 'l', "PROB", 0, "Lowest threshold of possibility of a data bin", 0},
    {"cut-x", 'x', "PROB", 0, "Cut the highest probability of met array", 0},
    {"cut-y", 'y', "PROB", 0, "Cut the highest probability of timing fluctuation array", 0},
    {"hz-ns", 'g', "FREQ", 0, "Tick per ns of the clock", 0},
    {"col-met", 'u', "COL", 0, "Valid data column in measurement csv file", 0},
    {"col-tf", 'v', "COL", 0, "Valid data column in timing fluctuation csv file", 0},
    {0}
};

// Parser function
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    filt_param_t *args = state->input;

    switch (key) {
        case 'm':
            args->in_tm_file = arg;
            break;
        case 's':
            args->in_tf_file = arg;
            break;
        case 'b':
            args->ref_time = atoi(arg);
            break;
        case 't':
            args->step = atoi(arg);
            break;
        case 'r':
            args->rnd = atoi(arg);
            break;
        case 'n':
            args->nsamp = atoi(arg);
            break;
        case 'l':
            args->p_low = atof(arg);
            break;
        case 'x':
            args->p_xcut = atof(arg);
            break;
        case 'y':
            args->p_ycut = atof(arg);
            break;
        case 'g':
            args->hz_ns = atof(arg);
            break;
        case 'u':
            args->met_col = atoi(arg);
            break;
        case 'v':
            args->tf_col = atoi(arg);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

// argp parser structure
static struct argp argp = {options, parse_opt, NULL, NULL};


/*=== BEGIN: Implementations ===*/

int
cmp(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}


int
sim_met(prob_hist_t *tmh, prob_hist_t *trh, i64 *tf_arr, u64 tf_len, double *psim, filt_param_t *args) {
    u64 nbin = tmh->nbin;
    double *ptrhc; // Cumulative probabilities of trh.
    struct timeval tv;

    ptrhc = (double *)malloc(nbin*sizeof(double));

    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1e6 + tv.tv_usec);

    // Calculating tr's cumulative probability
    ptrhc[0] = 0.0;
    psim[0] = 0;
    for (u64 ir = 1; ir < nbin; ir ++) {
        ptrhc[ir] = ptrhc[ir-1] + trh->pbin[ir-1].p;
        psim[ir] = 0;
    }
    
    for (u64 isamp = 0; isamp < args->nsamp; isamp ++){
        register double px = (double)rand() / (double)RAND_MAX;
        register u64 itf = rand() % tf_len;
        if (px < ptrhc[nbin-1]) {
            for (u64 ir = 0; ir < nbin - 1; ir ++) {
                if (px < ptrhc[ir+1]) {
                    register i64 tsim;
                    register double p0, p1, t0, t1, tx;
                    p0 = trh->pbin[ir].p;
                    p1 = trh->pbin[ir+1].p;
                    t0 = (double)trh->pbin[ir].t;
                    t1 = (double)trh->pbin[ir+1].t;
                    tx = t0 + (p1 - p0) / ((t1 - t0) * (px - p0));
                    tsim = (i64)tx + tf_arr[itf];
                    register u64 im = 0, ims = 0, ime = nbin - 2;
                    for (u64 im = 0; im < nbin - 1; im ++) {
                        if (tsim < tmh->pbin[im+1].t) {
                            psim[im] ++;
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }

    for (u64 isim = 0; isim < nbin; isim ++) {
        psim[isim] = psim[isim] / (double)args->nsamp;
    }

    free(ptrhc);
}

int 
calc_tr(prob_hist_t *tmh, prob_hist_t *tfh, prob_hist_t *trh, i64 *tf_arr, 
	u64 tf_len, filt_param_t *args) {
		
	u64 nbin, imb;
	i64 step;
	int err = 0;
	double tot_p = 0.0;	// Tracking total probability in trh for the exit condition
	double simh[tmh->nbin], simh2[tmh->nbin];	// Simulated probabilities

	// Init trh
	nbin = tmh->nbin;
	trh->nbin = nbin;
	trh->pbin = (prob_bin_t *)malloc(nbin * sizeof(prob_bin_t));
	if (trh->pbin == NULL) {
        printf("[FilT-calc_tr] ERROR, trh->pbin allocation failed.\n");
        err = errno;
        return err;
    }
	for (u64 i; i < nbin; i ++) {
		trh->pbin[i].t = tmh->pbin[i].t - tfh->pbin[0].t;
        // printf("%lld, %lld, %lld --- ", trh->pbin[i].t, tmh->pbin[i].t, tfh->pbin[0].t);
		trh->pbin[i].p = 0.0;
	}

	// Preparing auxiliaries.
	step = args->step;	// Optimization granularity in nanosecond.
    for (u64 i=0; i < nbin; i ++) {
        simh[i] = 0;
    }
	

	// Loop over each time section in trh
    imb = 0;
	while (imb < nbin - 1) {
		i64 tr_l = trh->pbin[imb].t, tr_r = trh->pbin[imb+1].t; // tr in [tr_l, tr_r)
        double p, tr_p = 0;
        double dp, dmdr, dp0, dp1;

        // ====== S1: Initial condition ======
        // Simulated prob > met in this bin, pass.
        p = tmh->pbin[imb].p - simh[imb];
        if (p <= 0) {
            trh->pbin[imb].p = 0;
            imb += 1;
            continue;
        }
        // Loop over each time step.
        for (i64 tr = tr_l; tr < tr_r; tr += step) {
            register i64 tf_l, tf_r;
            tf_l = tfh->pbin[0].t; 
            tf_r = tmh->pbin[imb+1].t-tr;
            for (u64 itf = 0; itf < tfh->nbin-1; itf ++) {
                if (tfh->pbin[itf+1].t > tf_r) {
                    //tr_p += (tf_r - tfh->pbin[itf].t) / (tfh->pbin[itf+1].t - tfh->pbin[itf].t) * tfh->pbin[itf].p;
                    break;
                }
                tr_p += tfh->pbin[itf].p;
                // printf("%f, %f, %f\n",p, tr_p, (double)step/(double)(tr_r - tr_l));
            }
        }
        if (tr_p < 1e-6) {
            tr_p = 0.001;
        } else {
            tr_p = p / tr_p / (double)((double)step / (double)(tr_r - tr_l));
        }
        trh->pbin[imb].p = tr_p;
        printf("[FilT-calc_tr] Init GUESS IMB=%lld, tr in [%lld, %lld), tr_p=%f\n", imb, tr_l, tr_r, tr_p);
        fflush(stdout);

        // ====== S2: Optimization ======
        sim_met(tmh, trh, tf_arr, tf_len, simh2, args);
        //for (u64 isim = 0; isim < nbin; isim ++) printf("%lld, %e\n", tmh->pbin[isim].t, simh2[isim]);
        dp = simh2[imb] - tmh->pbin[imb].p;
        dmdr = (simh2[imb] - simh[imb]) / tr_p;
        dp0 = 0x7fffffff;
        dp1 = dp;
        do {
            dp0 = dp1;
            tr_p = trh->pbin[imb].p;
            trh->pbin[imb].p -= dp1 / dmdr;
            // printf("[FilT-calc_tr] IMB=%lld, tr_p=%f, dp=%f\n", imb, tr_p, dp1);
            for (u64 isimh = 0; isimh < nbin; isimh ++) {
                simh[isimh] = simh2[isimh];
            }
            if (trh->pbin[imb].p <= 0) {
                // We do not want nagative probability
                dp1 = 0; // Don't sim, exit.
                break;
            } else {
                sim_met(tmh, trh, tf_arr, tf_len, simh2, args);
                dp1 = simh2[imb] - tmh->pbin[imb].p;
            }
        } while (fabs(dp0) > fabs(dp1));

        // exit condition
        if (tot_p + tr_p >= 1.0) {
            if (abs(tot_p + tr_p - 1) >= abs(tot_p - 1)) {
                tr_p = 0;
            } 
            trh->pbin[imb].p = tr_p;
            tot_p += tr_p;
            printf("[FilT-calc_tr] IMB=%lld, tr in [%lld, %llu), p=%f, tot_p=%f\n", 
                imb, trh->pbin[imb].t, trh->pbin[imb+1].t, trh->pbin[imb].p, tot_p);
            break;
        } else {
            trh->pbin[imb].p = tr_p;
            tot_p += tr_p;
            printf("[FilT-calc_tr] IMB=%lld, tr in [%lld, %llu), p=%f, tot_p=%f\n", 
                imb, trh->pbin[imb].t, trh->pbin[imb+1].t, trh->pbin[imb].p, tot_p);
            fflush(stdout);
            imb ++;
        }
	}
    printf("[FilT-calc_tr] Normalizing probabilities...");
    fflush(stdout);
    for (u64 i = 0; i < trh->nbin; i ++) {
        trh->pbin[i].p /= tot_p;
    //    printf("%lld, %e\n", trh->pbin[i].t, trh->pbin[i].p);
    } 
    printf("Done.\n");
    fflush(stdout);

	return err;
}

int
slice(i64 *arr, u64 len, i64 rnd, double p_low, double p_cut, prob_hist_t *phist) {
    struct timeval st, en;
    double duration, p;
    int err = 0;
    i64 tmin, tmax, t_bin_min, t_bin_max;
    u64 nbin, nc_low, nbv, nc; // init bin count, lowest point count, valid nbin
    u64 len2;    // Cut off the largest p_cut
    i64 *pbins;     // Left edge of each bin.
    u64 *pcnts, *pvcnts;    // Point counts of each bin.

    len2 = (u64)((double)len * (1.0 - p_cut)); // Cut off the largest p_cut
    // Ascending sort
    printf("[FilT-slice] Sorting ... ");
    gettimeofday(&st, NULL);
    qsort(arr, len, sizeof(i64), cmp);
    gettimeofday(&en, NULL);
    duration = (double)en.tv_sec + (double)(en.tv_usec * 1e-6) - \
               ((double)st.tv_sec + (double)(st.tv_usec * 1e-6));
    printf("%lld, %lld, ... %lld ", arr[0], arr[1], arr[len-1]);
    printf("done, in %.3fs.\n", duration);
    printf("[FilT-slice] len=%llu, len2=%llu\n", len, len2);

    // Counting
    tmin = arr[0];
    tmax = arr[len2-1];
    // Corner case, stable overhead.
    if (tmin == tmax) {
        t_bin_min = tmin;
        t_bin_max = tmax;
    } else {
        t_bin_min = tmin / rnd * rnd;
        t_bin_max = tmax / rnd * rnd + rnd;
    }
    nbin = (t_bin_max - t_bin_min) / rnd + 1;
    printf("[FilT-slice] nbin=%llu, tmin=%lld, t_bin_min=%lld, tmax=%lld, t_bin_max=%lld\n", 
        nbin, tmin, t_bin_min, tmax, t_bin_max);

    // Corner case, only 1 or 2 bin
    if (nbin <= 2) {
        phist->nbin = nbin;
        phist->pbin = (prob_bin_t *)malloc(nbin * sizeof(prob_bin_t));
        if (phist->pbin == NULL) {
            printf("[FilT-slice] ERROR, phist->pbin allocation failed.\n");
            err = errno;
            return err;
        }
		phist->pbin[nbin-1].t = t_bin_max;
        phist->pbin[nbin-1].p = 0;
        phist->pbin[0].t = t_bin_min;
        phist->pbin[0].p = 1.0;

		return err;
    }
    
    printf("[FilT-slice] Init bin count: %llu\n", nbin);
    pbins = (i64 *)malloc(nbin * sizeof(i64));
    pcnts = (u64 *)malloc(nbin * sizeof(u64));
    pvcnts = (u64 *)malloc(nbin * sizeof(u64));
    if (pbins == NULL || pcnts == NULL || pvcnts == NULL) {
        printf("[FilT-slice] ERROR, pbins and pcnts allocation failed.\n");
        err = errno;
        return err;
    }
    // Initalizing hist
    printf("[FilT-slice] Binning ... ");
    for (u64 i = 0; i < nbin; i ++) {
        pbins[i] = t_bin_min + i * rnd;
        pcnts[i] = 0;
    }

    // Counting in even bins
    for (u64 i = 1; i < len2; i ++) {
        u64 ib = (arr[i] - t_bin_min) / rnd;
        u64 ib_l = (arr[i-1] - t_bin_min) / rnd;
        if (ib_l != ib) {
            pcnts[ib_l] = i - 1;
        }
    }
    pcnts[(u64)((arr[len2-1]-t_bin_min)/rnd)] = len2; // Last bin

    // Merging for p_low
    nc_low = (u64)((double)len2 * p_low);
    nbv = 0;
    nc = 0;
	// The last bin is always empty, only for marking the right edge of the histo.
    for (u64 i = 0; i < nbin-1; i ++) {
        if (i == 0) {
            nc += pcnts[i];
        } else {
            nc += pcnts[i] - pcnts[i-1];
        }
        if (nc >= nc_low) {
            pbins[nbv+1] = pbins[i+1];
            pvcnts[nbv+1] = 0;
            pvcnts[nbv] = nc;
            nc = 0;
            nbv ++;
        }
    }
    pbins[nbv+1] = t_bin_max;
    pvcnts[nbv+1] = 0;
    pvcnts[nbv] = nc;
    nbv ++;

    // Calculating probability
    phist->nbin = nbv+1;
    phist->pbin = (prob_bin_t *)malloc((nbv + 1) * sizeof(prob_bin_t));
    if (phist->pbin == NULL) {
        printf("[FilT-slice] ERROR, phist->pbin allocation failed.\n");
        err = errno;
        return err;
    }
    for (u64 i = 0; i < nbv; i ++) {
        phist->pbin[i].t = pbins[i];
        phist->pbin[i].p = (double)pvcnts[i] / (double)len2;
    }
    phist->pbin[nbv].t = t_bin_max;
    phist->pbin[nbv].p = 0;
    printf("Done.\n");
    printf("[FilT-slice] Final bin count: %llu\n", phist->nbin);

    free(pbins);
    free(pcnts);
    free(pvcnts);

    return err;
}

int
parse_csv(char *fpath, int col, i64 adj, i64 rnd, double hz_ns, double p_low, double p_cut, 
    u64 *len, i64 **arr,  prob_hist_t *hist){

    FILE *fp;
    i64 fsize;
    i64 nrow;
    char *buf, *ptok;
    int err = 0;
    i64 i = 0;
    i64 icol = 0, irow = 0, irow_st = 0;

    printf("[FilT-parse_csv] Reading %s.\n", fpath);
    fp = fopen(fpath, "r");
    if (fp == NULL) {
        printf("[FilT-parse_csv] ERROR, file does not exist at given path: %s.\n", fpath);
        return errno;
    }
    fseek(fp, 0, SEEK_END);
    fsize = ftell(fp);
    rewind(fp);

    // Allocate memory for raw file.
    buf = (char *)malloc(fsize + 1);
    if (buf == NULL) {
        printf("[FilT-parse_csv] ERROR, read buffer allocation failed.\n");
        fclose(fp);
        return errno;
    }

    // Read the file content into the character array
    fread(buf, fsize, 1, fp);
    buf[fsize] = '\0';  // null terminator
    fclose(fp);

    
    // Counting how many column in a line.
    while (buf[i]!='\n') {
        if (buf[i] == ','){
            icol ++;
        }
        i ++;
    }
    // The data column exceeds existing column. error.
    if (col > icol) {
        printf("[FilT-parse_csv] ERROR, specified data column %d exceeds "
               "the valid column range [0, %lld].\n", col, icol);
        err = -1;
        return err;
    }

    // Counting how many rows in the file.
    i = 0;
    irow = 0;
    irow_st = 0;
    while (buf[i]!='\0') {
        if (buf[i] == '\n'){
            irow ++;
        }
        i ++;
    }
    nrow = irow;
    *len = (u64)((double)nrow * (1.0 - p_cut));
    printf("[FilT-parse_csv] %lld records are found, now extracting column %d.\n", nrow, col);

    // Allocating memory space for raw data array.
    *arr = (i64 *)malloc(nrow*sizeof(i64));
    if (*arr == NULL) {
        printf("[FilT-parse_csv] ERROR, read buffer allocation failed.\n");
        fclose(fp);
        err = errno;
        return err;
    }

    // Converting char to i64 at given column;
    
    i = 0;
    irow_st = 0;
    irow = 0;
    while (buf[i]!='\0') {
        // Loop over lines.
        if (buf[i] == '\n') {
            buf[i] = '\0';
            icol = 0;
            ptok = strtok(&buf[irow_st], ",");
            while (ptok != NULL) {
                if (icol == col) {
                    (*arr)[irow] = strtoll(ptok, NULL, 10);
                    irow ++;
                    break;
                }
                ptok = strtok(NULL, ",");
                icol ++;
            }
            irow_st = i + 1;
        }
        i ++;
    }
    // Converting tick to ns
    printf("[FilT-parse_csv] FREQ=%f, converting tick to ns.\n", hz_ns);
    for (u64 ir = 0; ir < nrow; ir ++) {
        fflush(stdout);
        register double tick = (double)(*arr)[ir];
        tick = tick / hz_ns + adj;
        (*arr)[ir] = (i64)tick;
        fflush(stdout);
    }

    // Slicing
    printf("[FilT-parse_csv] Start data binning.\n");
    err = slice(*arr, nrow, rnd, p_low, p_cut, hist);
    if (err != 0) {
        printf("[FilT-parse_csv] ERROR, slicing failed.\n");
        return err;
    }

    return 0;
}

int
run_filt(filt_param_t *args, prob_hist_t *tr_hist){
    u64 tm_len, tf_len;
    i64 *tm_arr, *tf_arr;
    prob_hist_t tm_hist, tf_hist;
    int err;

    // Parsing csv files and slicing specified column into histogram, saving to pmet_hist and ptf_hist
    // malloc inside slice function
    printf("[FilT-run_filt] Parsing measurement file %s\n", args->in_tm_file);
    err = parse_csv(args->in_tm_file, args->met_col, 0, args->rnd, args->hz_ns, 
                    args->p_low, args->p_xcut, &tm_len, &tm_arr, &tm_hist);
    if (err) {
        printf("[FilT-run_filt] Error in parsing measurement csv file. ERRCODE %d\n", err);
        return err;
    }
    //FILE *fp = fopen("met_hist.csv", "w");
    //for (i64 i = 0; i < tm_hist.nbin; i ++) {
    //    fprintf(fp, "%lld, %lf\n", tm_hist.pbin[i].t, tm_hist.pbin[i].p);
    //}
    //fclose(fp);

    printf("[FilT-run_filt] Parsing timing fluctuation file %s\n", args->in_tf_file);
    err = parse_csv(args->in_tf_file, args->tf_col, 0 - args->ref_time, args->rnd, args->hz_ns, 
                    args->p_low, args->p_ycut, &tf_len, &tf_arr, &tf_hist);
    if (err) {
        printf("[FilT-run_filt] Error in parsing timing fluctuations csv file. ERRCODE %d\n", err);
        return err;
    }

    //fp = fopen("tf_hist.csv", "w");
    //for (i64 i = 0; i < tf_hist.nbin; i ++) {
    //    fprintf(fp, "%lld, %lf\n", tf_hist.pbin[i].t, tf_hist.pbin[i].p);
    //}
    //fclose(fp);

	if (tf_hist.nbin <= 2) {
		printf("[FilT-run_filt] Warning: few timing fluctuation deceted, will not run FilT.\n");
	} else {
        printf("[FilT-run_filt] Start estimating real run time distribution.\n");
        err = calc_tr(&tm_hist, &tf_hist, tr_hist, tf_arr, tf_len, args);
        if (err) {
            printf("[FilT-run_filt] Error in transposed convolution. ERRCODE %d\n", err);
            return err;
        }
	}

    free(tm_arr);
    free(tf_arr);
    free(tm_hist.pbin);
    free(tf_hist.pbin);
    return 0;
}

//int
//arg_parse(int argc, char **argv){
//
//}

/*=== END: Implementations ===*/

/*=== BEGIN: Main Entry ===*/

int
main(int argc, char **argv){
    filt_param_t args;
    prob_hist_t tr_hist;
    i64 res_len;
    int err;
    
    args.in_tm_file = "met.csv";
    args.in_tf_file = "tf.csv";
    args.out_tr_file = "filt_out.csv";
    args.ref_time = 90000;
    args.rnd = 100;
    args.nsamp = 1000000;
    args.nsim = 1;
    args.step = 10;
    args.met_col = 2;
    args.tf_col = 2;
    args.p_low = 0.001;
    args.p_xcut = 0.001;
    args.p_ycut = 0.001;
    args.hz_ns = 1.0;

    // Parse command line
    argp_parse(&argp, argc, argv, 0, 0, &args);

    // Reading input files and estimating the real run time distribution.
    err = run_filt(&args, &tr_hist);

    if (err == 0) {
        printf("[FilT-main] run_filt has been Finished. Writing results to %s\n", args.out_tr_file);
        FILE *fp = fopen(args.out_tr_file, "w");
        for (i64 i = 0; i < tr_hist.nbin; i ++) {
            fprintf(fp, "%lld, %lf\n", tr_hist.pbin[i].t, tr_hist.pbin[i].p);
        }
        fclose(fp);
        free(tr_hist.pbin);
    } else {
        printf("[FilT-main] run_filt returned with errors. ERRCODE %d\n", err);
        return err;
    }

    return 0;
}

/*=== END: Main Entry ===*/
