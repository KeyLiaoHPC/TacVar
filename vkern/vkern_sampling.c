#define _GNU_SOURCE
#define _ISOC11_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sched.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "mpi.h"

#ifdef USE_PAPI
#include "papi.h"
#endif

// Warmup for 1000ms.
#ifndef NWARM
#define NWARM 1000
#endif

// Ignoring some tests at the beginning
#ifndef NPASS
#define NPASS 2
#endif

// Number of tests for each interval
#ifndef NTEST
#define NTEST 10
#endif

#ifndef TBASE
#define TBASE 100
#endif

// Uniform: V1: Number of DSub between two bins; V2: Number of bins.
// Normal: V1: sigma, standard error
// Pareto: V1: alpha, Pareto exponent 
#ifndef V1
#define V1 100
#endif

// Number of intervals (>=1)
#ifndef V2
#define V2 20
#endif

#ifndef FSIZE
#define FSIZE 0
#endif

#ifndef LCUT
#define LCUT 1-0x7fffffff
#endif

#ifndef HCUT
#define HCUT 0x7fffffff-1
#endif

// Prefetching the benchmarking instructions.
#ifndef NPRECALC
#define NPRECALC 1
#endif

#ifndef NSAMP
#define NSAMP 0
#endif

// Timing macros
#define _read_ns(_ns) \
    do {                                                \
        register uint64_t ns;                           \
        asm volatile(                                   \
            "\n\tRDTSCP"                                \
            "\n\tshl $32, %%rdx"                        \
            "\n\tor  %%rax, %%rdx"                      \
            "\n\tmov %%rdx, %0"                         \
            "\n\t"                                      \
            :"=r" (ns)                                  \
            :                                           \
            : "memory", "%rax", "%rdx");                \
        _ns = ns;                                       \
    } while(0);

#define _mfence asm volatile("lfence"   "\n\t":::);


/**
 * @brief A desginated one-line computing kernel.
 */
void
flush_cache(double *pf_a, double *pf_b, double *pf_c, uint64_t npf) {
#ifdef INIT
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = i;
    }

#elif TRIAD
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = 0.42 * pf_b[i] + pf_c[i];
    }

#elif SCALE
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = 1.0001 * pf_b[i];
        pf_b[i] = 1.0001 * pf_a[i];
    }

#endif
    return;
}

/**
 * @brief Reading urandom file for random seeds.
 */
int 
get_urandom(uint64_t *x) {
    FILE *fp = fopen("/dev/urandom", "r");
    if (fp == NULL) {
        printf("Failed to open random file.\n");
        return 1;
    } 

    fread(x, sizeof(uint64_t), 1, fp); 

    fclose(fp);

    return 0;
}

/**
 * @brief Generate a shuffled list to loop through all ADD array's lengths.
 */
void
gen_walklist(uint64_t *len_list) {
#ifndef WALK_FILE
#ifdef UNIFORM
    for (int it = 0; it < NTEST; it ++) {
        for (int i = 0; i < V2; i ++) {
            len_list[NPASS+it*V2+i] = TBASE + i * V1;
        }
    }
    // Shuffle
    // Random seeds using nanosec timestamp
    struct timespec tv;
    uint64_t sec, nsec;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    sec = tv.tv_sec;
    nsec = tv.tv_nsec;
    nsec = sec * 1e9 + nsec + NWARM * 1e6;
    srand(nsec); 
    for (int i = NPASS; i < NPASS + NTEST * V2; i ++){
        int r = (int) (((float)rand() / (float)RAND_MAX) * (float)(NTEST * V2));
        uint64_t temp = len_list[i];
        // Randomly swapping two vkern loop length
        len_list[i] = len_list[NPASS+r];
        len_list[NPASS+r] = temp;
    }
    for (int i = 0; i < NPASS; i ++) {
        len_list[i] = TBASE + V1 * V2;
    }
#else
    double v1 = V1;
    uint64_t seed;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlux389;
    r = gsl_rng_alloc(T);
    get_urandom(&seed);
    gsl_rng_set(r, seed);
    for (int i = 0; i < NPASS+NTEST; i++) {
        double x = 0;
#ifdef NORMAL
        do {
            x = 1.0 + gsl_ran_gaussian(r, v1);
        } while (x < LCUT || x > HCUT);
#elif PARETO
        do {
            x = gsl_ran_pareto(r, v1, 1.0);
        } while (x > HCUT);
#endif
        len_list[i] = (int)(x*TBASE);
    }
    gsl_rng_free(r);
#endif
#else
    // Reading walk list from file
    FILE *fp = fopen("walk.csv", "r");
    for (int it = 0; it < NTEST + NPASS; it ++) {
        fscanf(fp, "%llu\n", len_list+it);
    }
    fclose(fp);
#endif
}

int
main(int argc, char **argv) {
    int ntest;
    register uint64_t a, b, c, d;
    uint64_t *p_len, *p_cy, *p_ns;
    uint64_t register cy0, cy1, ns0, ns1;
    uint64_t mask = (1 << 30), tbase=TBASE, fsize = FSIZE, npf = 0; // npf: flush arr length
    double *pf_a, *pf_b, *pf_c;
    int myrank, mycpu, nrank, errid;
#ifdef UNIFORM
    uint64_t v1 = V1, v2 = V2;

#elif NORMAL
    double v1 = V1;

#elif PARETO
    double v1 = V1;
#endif

    errid = MPI_Init(NULL, NULL);
    if (errid != MPI_SUCCESS) {
        printf("Faild to init MPI.\n");
        exit(1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    mycpu = sched_getcpu();

#ifdef UNIFORM
    ntest = NPASS + NTEST * V2;
#else
    ntest = NTEST + NPASS;
#endif

#ifdef USE_PAPI
    // Init PAPI
    int eventset = PAPI_NULL;
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&eventset);
    PAPI_start(eventset);
#endif

    p_len = (uint64_t *)malloc(ntest * sizeof(uint64_t));
    p_ns = (uint64_t *)malloc(ntest * sizeof(uint64_t));
    if (fsize) {
        npf = FSIZE / sizeof(double);
        pf_a = (double *)malloc(npf * sizeof(double));
        pf_b = (double *)malloc(npf * sizeof(double));
        pf_c = (double *)malloc(npf * sizeof(double));
    }
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = 1.1;
        pf_b[i] = 1.1;
        pf_c[i] = 1.1;
    }

    if (myrank == 0) {
#ifdef UNIFORM
        printf("Generating uniform ditribution. Tbase=%llu, interval=%llu, nint=%llu, Ntest=%llu.\n",
            tbase, v1, v2, NTEST);

#elif NORMAL
        printf("Generating normal ditribution. Tbase=%llu, sigma=%.4f, Ntest=%llu.\n",
            tbase, v1, NTEST);

#elif PARETO
        printf("Generating Pareto ditribution. Tbase=%llu, alpha=%.4f, Ntest=%llu.\n",
            tbase, v1, NTEST);
#else
        printf("Reading vkern walk list from walk.csv.\n");

#endif
        fflush(stdout);
        gen_walklist(p_len);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(p_len, ntest, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    a = 0;
    b = 0;
    if (argc < 10) {
        c = 1;
    } else {
        c = atoi(argv[1]);
    }
    d = 0;

    // Warm up
    do {
        if (myrank == 0) {
            printf("Warming up for %d ms.\n", NWARM);
        }
        struct timespec tv;
        uint64_t volatile sec, nsec; // For warmup
        clock_gettime(CLOCK_MONOTONIC, &tv);
        sec = tv.tv_sec;
        nsec = tv.tv_nsec;
        nsec = sec * 1e9 + nsec + NWARM * 1e6;

        while (tv.tv_sec * 1e9 + tv.tv_nsec < nsec) {
            a = b + 1;
            b = a + 1;
            clock_gettime(CLOCK_MONOTONIC, &tv);
        }
    } while (0);

    if (myrank == 0) {
        printf("Start random walking.\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int iwalk = 0; iwalk < ntest; iwalk ++) {
        register uint64_t n = p_len[iwalk] * 2;
        register uint64_t npre = NPRECALC;
        register uint64_t ra, rb, rc;
        struct timespec tv;

        // Flushing
        if (npf) {
            flush_cache(pf_a, pf_b, pf_c, npf);
        }

        // Instruction preload.
        ra = n + npre * 2;
        rc = c;
        while (ra != n) {
            rb = ra - rc;
            ra = rb - rc;
        }

        // Timing.
        while (ra!= 0) {
            rb = ra - rc;
            ra = rb - rc;
        }

        ra = NSAMP * 2;
#ifdef USE_PAPI
        ns0 = PAPI_get_real_nsec();

#elif USE_CGT
        clock_gettime(CLOCK_MONOTONIC, &tv);
        ns0 = tv.tv_sec * 1e9 + tv.tv_nsec;
#elif USE_WTIME
        ns0 = (uint64_t)(MPI_Wtime() * 1e9);
#else
        _read_ns (ns0);
        _mfence;
#endif
        while (ra!= 0) {
            rb = ra - rc;
            ra = rb - rc;
        }
#ifdef USE_PAPI
        ns1 = PAPI_get_real_nsec();

#elif USE_CGT
        clock_gettime(CLOCK_MONOTONIC, &tv);
        ns1 = tv.tv_sec * 1e9 + tv.tv_nsec;
#elif USE_WTIME
        ns1 = (uint64_t)(MPI_Wtime() * 1e9);
#else
        _read_ns (ns1);
        _mfence;
#endif

        p_ns[iwalk] = ns1 - ns0;
        MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef USE_PAPI
    PAPI_shutdown();
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    uint64_t *p_cy_all, *p_ns_all;
    p_ns_all = (uint64_t *)malloc(nrank*ntest*sizeof(uint64_t));
    MPI_Gather(p_ns, ntest, MPI_UINT64_T, p_ns_all, ntest, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
#ifdef USE_PAPI
        FILE *fp = fopen("papi_time_sample.csv", "w");
#elif USE_CGT
        FILE *fp = fopen("cgt_time_sample.csv", "w");
#elif USE_WTIME
        FILE *fp = fopen("wtime_time_sample.csv", "w");
#else
        FILE *fp = fopen("stiming_time_sample.csv", "w");
#endif
        for (int irank = 0; irank < nrank; irank ++) {
            for (int iwalk = NPASS; iwalk < ntest; iwalk ++) {
                int iall = irank * ntest + iwalk;
                fprintf(fp, "%d,%llu,%llu\n", irank, NSAMP, p_ns_all[iall]);
            }
        }
        fclose(fp);
    }
    free(p_ns_all);

    MPI_Finalize();

    if (npf) {
        free(pf_a);
        free(pf_b);
        free(pf_c);
    }

    free(p_len);
    free(p_ns);

    return 0;
}
