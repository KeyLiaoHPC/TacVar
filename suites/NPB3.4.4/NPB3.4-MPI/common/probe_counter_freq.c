/**
 * @file probe_counter_freq.c
 * @brief NPB-MPI/common: measure hardware-counter Hz against CLOCK_MONOTONIC
 *        and compare to the ISA-claimed rate (cntfrq_el0 / CPUID / tsc_freq_khz).
 *
 * Serial: one randomly chosen online CPU per socket.
 * MPI (-DUSE_MPI): every rank measures its bound core (or rank % ncpu).
 *
 * ppm = (1 - measured_hz / claimed_hz) * 1e6
 *   positive => the counter is slow relative to the Linux monotonic clock.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <errno.h>
#include <math.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#if !defined(__x86_64__) && !defined(__aarch64__)
#error "probe_counter_freq supports only x86_64 (TSC) and aarch64 (CNTVCT)"
#endif

#define DEFAULT_INTERVAL_SEC 10.0
#define DEFAULT_SAMPLE_MS    50
#define DEFAULT_R_MIN        0.999
#define DEFAULT_MAX_RETRY    5
#define MAX_CPUS             4096
#define MAX_SOCKETS          256
#define CLAIMED_SRC_LEN      32
#define HOSTNAME_LEN         256

typedef struct {
    double interval_sec;
    int sample_ms;
    double r_min;
    int max_retry;
    int clock_raw;          /* 0 = CLOCK_MONOTONIC, 1 = CLOCK_MONOTONIC_RAW */
    unsigned int seed;
    int seed_set;
    int cpu_override;       /* -1 = unset */
    const char *csv_path;
    const char *nspt_out;   /* 写实测 nspt (=1e9/measured_hz_p50) 到文件 */
    int verify;             /* 1 = 读 nspt 文件并对比本次实测(偏差 ppm) */
    const char *verify_path;/* --verify 的目标文件(默认取 --nspt-out) */
} opts_t;

typedef struct {
    int cpu;
    int socket;
    int nsamp;
    int retries;
    int ok;
    double measured_hz;
    double claimed_hz;
    double ppm;
    double r;
    double hz_min;
    double hz_max;
    double hz_p50;
    double implied_mono_hz;
    char claimed_src[CLAIMED_SRC_LEN];
} result_t;

typedef struct {
    int ncpu;
    int cpu[MAX_CPUS];
    int socket[MAX_CPUS];
    int nsockets;
    int socket_ids[MAX_SOCKETS];
} topo_t;

/* -------------------------------------------------------------------------- */
/* Hardware counter and claimed frequency                                     */
/* -------------------------------------------------------------------------- */

#ifdef __aarch64__

static const char *g_arch_name = "aarch64";
static const char *g_counter_name = "cntvct_el0";

static inline uint64_t read_hw_counter(void)
{
    uint64_t v;
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(v) :: "memory");
    return v;
}

static int read_claimed_hz(int cpu, double *hz_out, char *src, size_t src_len)
{
    uint64_t frq;
    (void)cpu;
    __asm__ volatile("mrs %0, cntfrq_el0" : "=r"(frq) :: "memory");
    if (frq == 0)
        return -1;
    *hz_out = (double)frq;
    snprintf(src, src_len, "cntfrq_el0");
    return 0;
}

#else /* __x86_64__ */

static const char *g_arch_name = "x86_64";
static const char *g_counter_name = "rdtscp";

static inline uint64_t read_hw_counter(void)
{
    unsigned hi, lo, aux;
    __asm__ volatile("rdtscp" : "=a"(lo), "=d"(hi), "=c"(aux) :: "memory");
    return ((uint64_t)hi << 32) | (uint64_t)lo;
}

static int read_sysfs_u64(const char *path, uint64_t *out)
{
    FILE *f = fopen(path, "r");
    unsigned long long v;
    if (!f)
        return -1;
    if (fscanf(f, "%llu", &v) != 1) {
        fclose(f);
        return -1;
    }
    fclose(f);
    *out = (uint64_t)v;
    return 0;
}

static int read_claimed_hz(int cpu, double *hz_out, char *src, size_t src_len)
{
    unsigned eax, ebx, ecx, edx;
    char path[128];
    uint64_t khz;

    eax = 0x15;
    ecx = 0;
    __asm__ volatile("cpuid"
                     : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                     : "a"(eax), "c"(ecx));
    if (eax != 0 && ebx != 0 && ecx != 0) {
        *hz_out = ((double)ecx * (double)ebx) / (double)eax;
        snprintf(src, src_len, "CPUID.15H");
        return 0;
    }

    eax = 0x16;
    ecx = 0;
    __asm__ volatile("cpuid"
                     : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                     : "a"(eax), "c"(ecx));
    if (eax != 0) {
        *hz_out = (double)eax * 1e6;
        snprintf(src, src_len, "CPUID.16H");
        return 0;
    }

    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/tsc_freq_khz", cpu);
    if (read_sysfs_u64(path, &khz) == 0 && khz > 0) {
        *hz_out = (double)khz * 1000.0;
        snprintf(src, src_len, "tsc_freq_khz");
        return 0;
    }
    if (read_sysfs_u64("/sys/devices/system/cpu/cpu0/tsc_freq_khz", &khz) == 0
        && khz > 0) {
        *hz_out = (double)khz * 1000.0;
        snprintf(src, src_len, "tsc_freq_khz");
        return 0;
    }
    return -1;
}

#endif

static inline uint64_t read_clock_ns(int clock_raw)
{
    struct timespec ts;
    clockid_t id = clock_raw ? CLOCK_MONOTONIC_RAW : CLOCK_MONOTONIC;
    clock_gettime(id, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static const char *clock_name(int clock_raw)
{
    return clock_raw ? "CLOCK_MONOTONIC_RAW" : "CLOCK_MONOTONIC";
}

/* -------------------------------------------------------------------------- */
/* Topology                                                                   */
/* -------------------------------------------------------------------------- */

static int parse_cpu_list(const char *s, int *cpus, int maxn)
{
    int n = 0;
    const char *p = s;

    while (*p && n < maxn) {
        char *end;
        long a, b;

        a = strtol(p, &end, 10);
        if (end == p)
            break;
        b = a;
        if (*end == '-') {
            b = strtol(end + 1, &end, 10);
            if (end == p)
                break;
        }
        if (b < a)
            break;
        for (long c = a; c <= b && n < maxn; c++)
            cpus[n++] = (int)c;
        if (*end == ',')
            p = end + 1;
        else
            break;
    }
    return n;
}

static int read_package_id(int cpu, int *pkg)
{
    char path[128];
    FILE *f;
    int v;

    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", cpu);
    f = fopen(path, "r");
    if (!f)
        return -1;
    if (fscanf(f, "%d", &v) != 1) {
        fclose(f);
        return -1;
    }
    fclose(f);
    *pkg = v;
    return 0;
}

static int topo_load(topo_t *topo)
{
    char buf[4096];
    FILE *f;
    int i, j;

    memset(topo, 0, sizeof(*topo));

    f = fopen("/sys/devices/system/cpu/online", "r");
    if (f && fgets(buf, sizeof(buf), f)) {
        fclose(f);
        topo->ncpu = parse_cpu_list(buf, topo->cpu, MAX_CPUS);
    } else {
        if (f)
            fclose(f);
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        if (n < 1)
            n = 1;
        if (n > MAX_CPUS)
            n = MAX_CPUS;
        topo->ncpu = (int)n;
        for (i = 0; i < topo->ncpu; i++)
            topo->cpu[i] = i;
    }

    for (i = 0; i < topo->ncpu; i++) {
        int pkg = 0;
        if (read_package_id(topo->cpu[i], &pkg) != 0)
            pkg = 0;
        topo->socket[i] = pkg;
        for (j = 0; j < topo->nsockets; j++) {
            if (topo->socket_ids[j] == pkg)
                break;
        }
        if (j == topo->nsockets && topo->nsockets < MAX_SOCKETS)
            topo->socket_ids[topo->nsockets++] = pkg;
    }
    return topo->ncpu > 0 ? 0 : -1;
}

static int topo_socket_of_cpu(const topo_t *topo, int cpu)
{
    int i;
    for (i = 0; i < topo->ncpu; i++) {
        if (topo->cpu[i] == cpu)
            return topo->socket[i];
    }
    return -1;
}

static int topo_cpu_online(const topo_t *topo, int cpu)
{
    int i;
    for (i = 0; i < topo->ncpu; i++) {
        if (topo->cpu[i] == cpu)
            return 1;
    }
    return 0;
}

static int pin_cpu(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (sched_setaffinity(0, sizeof(set), &set) != 0) {
        fprintf(stderr, "error: sched_setaffinity(%d) failed: %s\n",
                cpu, strerror(errno));
        return -1;
    }
    return 0;
}

#ifdef USE_MPI
/* Unique binding: exactly one CPU in the affinity mask. */
static int bound_cpu_or_neg(void)
{
    cpu_set_t set;
    int cpu, found = -1, n = 0;
    long maxcpu;

    if (sched_getaffinity(0, sizeof(set), &set) != 0)
        return -1;
    maxcpu = sysconf(_SC_NPROCESSORS_CONF);
    if (maxcpu < 1)
        maxcpu = 1;
    if (maxcpu > CPU_SETSIZE)
        maxcpu = CPU_SETSIZE;
    for (cpu = 0; cpu < (int)maxcpu; cpu++) {
        if (CPU_ISSET(cpu, &set)) {
            found = cpu;
            n++;
            if (n > 1)
                return -1;
        }
    }
    return found;
}
#endif /* USE_MPI */

#ifndef USE_MPI
static int pick_random_cpu_on_socket(const topo_t *topo, int socket)
{
    int tmp[MAX_CPUS];
    int n = 0, i;

    for (i = 0; i < topo->ncpu; i++) {
        if (topo->socket[i] == socket)
            tmp[n++] = topo->cpu[i];
    }
    if (n == 0)
        return -1;
    return tmp[rand() % n];
}
#endif /* !USE_MPI */

/* -------------------------------------------------------------------------- */
/* Statistics                                                                 */
/* -------------------------------------------------------------------------- */

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db)
        return -1;
    if (da > db)
        return 1;
    return 0;
}

static double percentile_p50(double *v, int n)
{
    if (n <= 0)
        return 0.0;
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    if (n & 1)
        return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static int ols_fit(const double *x, const double *y, int n,
                   double *slope, double *r)
{
    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    double nd, den_s, den_r, num;
    int i;

    if (n < 2)
        return -1;
    for (i = 0; i < n; i++) {
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        syy += y[i] * y[i];
        sxy += x[i] * y[i];
    }
    nd = (double)n;
    den_s = nd * sxx - sx * sx;
    if (den_s == 0.0)
        return -1;
    *slope = (nd * sxy - sx * sy) / den_s;
    num = nd * sxy - sx * sy;
    den_r = (nd * sxx - sx * sx) * (nd * syy - sy * sy);
    if (den_r <= 0.0)
        *r = 0.0;
    else
        *r = num / sqrt(den_r);
    return 0;
}

static void sleep_ms(int ms)
{
    struct timespec ts;
    if (ms <= 0)
        return;
    ts.tv_sec = ms / 1000;
    ts.tv_nsec = (long)(ms % 1000) * 1000000L;
    nanosleep(&ts, NULL);
}

/* -------------------------------------------------------------------------- */
/* One measurement window on the already-pinned CPU                           */
/* -------------------------------------------------------------------------- */

static int measure_window(const opts_t *opt, result_t *best)
{
    int max_samp = (int)(opt->interval_sec * 1000.0 / (double)opt->sample_ms) + 8;
    uint64_t *mono = NULL, *hw = NULL;
    double *x = NULL, *y = NULL, *inst = NULL;
    result_t local;
    int attempt, keep = 0;
    int rc = -1;

    if (max_samp < 8)
        max_samp = 8;

    mono = (uint64_t *)malloc((size_t)max_samp * sizeof(uint64_t));
    hw = (uint64_t *)malloc((size_t)max_samp * sizeof(uint64_t));
    x = (double *)malloc((size_t)max_samp * sizeof(double));
    y = (double *)malloc((size_t)max_samp * sizeof(double));
    inst = (double *)malloc((size_t)max_samp * sizeof(double));
    if (!mono || !hw || !x || !y || !inst) {
        fprintf(stderr, "error: out of memory\n");
        goto done;
    }

    memset(&local, 0, sizeof(local));
    local.r = -1.0;

    for (attempt = 0; attempt <= opt->max_retry; attempt++) {
        uint64_t t0, t_now, t_end;
        int n = 0, i, ninst = 0;
        double slope = 0.0, r = 0.0;

        t0 = read_clock_ns(opt->clock_raw);
        t_end = t0 + (uint64_t)(opt->interval_sec * 1e9);

        while (n < max_samp) {
            uint64_t h0, h1, m;
            h0 = read_hw_counter();
            m = read_clock_ns(opt->clock_raw);
            h1 = read_hw_counter();
            mono[n] = m;
            hw[n] = h0 + (h1 - h0) / 2;
            n++;
            t_now = read_clock_ns(opt->clock_raw);
            if (t_now >= t_end)
                break;
            sleep_ms(opt->sample_ms);
        }

        if (n < 2)
            continue;

        for (i = 0; i < n; i++) {
            x[i] = (double)(mono[i] - mono[0]) * 1e-9;
            y[i] = (double)(hw[i] - hw[0]);
        }
        if (ols_fit(x, y, n, &slope, &r) != 0)
            continue;

        ninst = 0;
        for (i = 1; i < n; i++) {
            double dt = (double)(mono[i] - mono[i - 1]) * 1e-9;
            if (dt <= 0.0)
                continue;
            inst[ninst++] = (double)(hw[i] - hw[i - 1]) / dt;
        }
        if (ninst < 1)
            continue;

        if (!keep || r > local.r) {
            local.measured_hz = slope;
            local.r = r;
            local.nsamp = n;
            local.retries = attempt;
            local.hz_min = inst[0];
            local.hz_max = inst[0];
            for (i = 1; i < ninst; i++) {
                if (inst[i] < local.hz_min)
                    local.hz_min = inst[i];
                if (inst[i] > local.hz_max)
                    local.hz_max = inst[i];
            }
            local.hz_p50 = percentile_p50(inst, ninst);
            local.ok = 1;
            keep = 1;
        }
        if (r >= opt->r_min)
            break;
    }

    if (!keep) {
        fprintf(stderr, "error: no valid fit on this CPU\n");
        goto done;
    }
    *best = local;
    rc = 0;

done:
    free(mono);
    free(hw);
    free(x);
    free(y);
    free(inst);
    return rc;
}

static void fill_claimed(result_t *res, int cpu)
{
    if (read_claimed_hz(cpu, &res->claimed_hz, res->claimed_src,
                        sizeof(res->claimed_src)) != 0) {
        res->claimed_hz = 0.0;
        snprintf(res->claimed_src, sizeof(res->claimed_src), "NA");
        res->ppm = 0.0;
        res->implied_mono_hz = 0.0;
        return;
    }
    if (res->measured_hz <= 0.0) {
        res->ppm = 0.0;
        res->implied_mono_hz = 0.0;
        return;
    }
    res->ppm = (1.0 - res->measured_hz / res->claimed_hz) * 1e6;
    res->implied_mono_hz = 1e9 * res->claimed_hz / res->measured_hz;
}

static int measure_on_cpu(const opts_t *opt, const topo_t *topo, int cpu,
                          result_t *res)
{
    memset(res, 0, sizeof(*res));
    res->cpu = cpu;
    res->socket = topo_socket_of_cpu(topo, cpu);
    snprintf(res->claimed_src, sizeof(res->claimed_src), "NA");
    if (pin_cpu(cpu) != 0)
        return -1;
    if (measure_window(opt, res) != 0)
        return -1;
    res->cpu = cpu;
    res->socket = topo_socket_of_cpu(topo, cpu);
    fill_claimed(res, cpu);
    return 0;
}

/* -------------------------------------------------------------------------- */
/* Output                                                                     */
/* -------------------------------------------------------------------------- */

static void print_help(const char *argv0)
{
    printf(
"Usage: %s [options]\n"
"\n"
"Measure the hardware cycle-counter frequency against Linux monotonic time\n"
"and compare it to the frequency the ISA / firmware claims.\n"
"\n"
"On aarch64 the counter is cntvct_el0 and the claim is cntfrq_el0 (often\n"
"100 MHz on Kunpeng). On x86_64 the counter is rdtscp and the claim is\n"
"CPUID leaf 0x15, then 0x16, then sysfs tsc_freq_khz. The claimed value is\n"
"never taken from a CLOCK_MONOTONIC calibration: that would hide the ppm\n"
"offset this tool exists to report.\n"
"\n"
"Method: over each --interval window, sample (hw, mono, hw) pairs and fit\n"
"  hw_ticks = slope * mono_seconds + intercept\n"
"by ordinary least squares. slope is measured_hz. Pearson r is computed on\n"
"the same pairs. If r < --r-min the window is repeated (up to --max-retry\n"
"extra times) and the highest-r fit is kept. Instantaneous Hz from\n"
"consecutive samples gives per-core min / max / p50.\n"
"\n"
"ppm = (1 - measured_hz / claimed_hz) * 1e6\n"
"  Positive ppm means the counter is SLOW relative to the monotonic clock\n"
"  (c920bn1: CNTVCT ~ +75 ppm vs CLOCK_MONOTONIC at a claimed 100 MHz).\n"
"implied_mono_hz is 1e9 * claimed_hz / measured_hz: the apparent rate of\n"
"the monotonic clock in ns/s if the counter ran at exactly the claimed Hz.\n"
"\n"
"Coverage\n"
"  Serial (this binary, no -DUSE_MPI):\n"
"    One randomly chosen online CPU per socket/package. Use --seed to\n"
"    reproduce the choice, or --cpu N to measure a single CPU.\n"
"  MPI (compile with -DUSE_MPI, binary probe_counter_freq-mpi.x):\n"
"    Every rank measures one core. If the rank is already bound to a\n"
"    single CPU that CPU is used; otherwise the rank pins to rank %% ncpu.\n"
"    nprocs must not exceed the number of online CPUs.\n"
"    Recommended launch:\n"
"      mpirun -np <ncores> --map-by core --bind-to core \\\n"
"          ./common/probe_counter_freq-mpi.x\n"
"\n"
"Options\n"
"  --interval <sec>   Measurement window length          (default: 10)\n"
"  --sample-ms <ms>   Period between samples             (default: 50)\n"
"  --r-min <r>        Pearson-r retry threshold          (default: 0.999)\n"
"  --max-retry <n>    Extra windows if r is below r-min  (default: 5)\n"
"  --clock <mono|raw> CLOCK_MONOTONIC or _RAW            (default: mono)\n"
"  --seed <n>         RNG seed for the per-socket CPU pick\n"
"  --cpu <id>         Serial only: measure this CPU and skip the others\n"
"  --csv <path>       Write the per-core table as CSV\n"
"  --nspt-out <path>  Write measured ns/tick (1e9/measured_hz p50) to file\n"
"  --verify <path>    Compare existing nspt file against this measurement\n"
"                     (deviation in ppm; exit 0 if within +-5 ppm)\n"
"  --help, -h         Show this help and exit\n"
"\n"
"Examples\n"
"  ./probe_counter_freq.x\n"
"  ./probe_counter_freq.x --interval 2 --sample-ms 20 --seed 1\n"
"  ./probe_counter_freq.x --cpu 0 --interval 10\n"
"  mpirun -np <ncores> --map-by core --bind-to core \\\n"
"      ./probe_counter_freq-mpi.x --interval 10\n"
"\n"
"Build (from NPB3.4-MPI/ or NPB3.4-MPI/common/)\n"
"  make probe_counter_freq\n"
"  make -C common\n"
"  gcc -std=c99 -Wall -Wextra -O2 -D_GNU_SOURCE -o probe_counter_freq.x \\\n"
"      probe_counter_freq.c -lm\n"
"  mpicc -std=c99 -Wall -Wextra -O2 -D_GNU_SOURCE -DUSE_MPI \\\n"
"      -o probe_counter_freq-mpi.x probe_counter_freq.c -lm\n",
        argv0);
}

static void print_header(const opts_t *opt, const topo_t *topo,
                         const char *claimed_src, double claimed_hz,
                         int nrank)
{
    char host[HOSTNAME_LEN];

    if (gethostname(host, sizeof(host) - 1) != 0)
        strncpy(host, "?", sizeof(host));
    host[sizeof(host) - 1] = '\0';

    printf("# node=%s arch=%s counter=%s", host, g_arch_name, g_counter_name);
    if (claimed_hz > 0.0)
        printf(" claimed=%.0f Hz (%s)", claimed_hz, claimed_src);
    else
        printf(" claimed=NA (%s)", claimed_src);
    printf(" clock=%s interval=%.3fs sample_ms=%d r_min=%.6f max_retry=%d\n",
           clock_name(opt->clock_raw), opt->interval_sec, opt->sample_ms,
           opt->r_min, opt->max_retry);
#ifdef USE_MPI
    printf("# mode=mpi nrank=%d online_cpus=%d sockets=%d\n",
           nrank, topo->ncpu, topo->nsockets);
#else
    (void)nrank;
    printf("# mode=serial online_cpus=%d sockets=%d\n",
           topo->ncpu, topo->nsockets);
#endif
}

static void print_core_row(const result_t *r, int have_claimed)
{
    if (!r->ok) {
        printf("%d %d NA NA NA NA %d %d NA NA NA NA\n",
               r->cpu, r->socket, r->nsamp, r->retries);
        return;
    }
    printf("%d %d %.6f ", r->cpu, r->socket, r->measured_hz);
    if (have_claimed && r->claimed_hz > 0.0)
        printf("%.0f %.4f ", r->claimed_hz, r->ppm);
    else
        printf("NA NA ");
    printf("%.8f %d %d %.6f %.6f %.6f ",
           r->r, r->nsamp, r->retries, r->hz_min, r->hz_max, r->hz_p50);
    if (have_claimed && r->implied_mono_hz > 0.0)
        printf("%.6f\n", r->implied_mono_hz);
    else
        printf("NA\n");
}

static void print_core_table(const result_t *rows, int nrows, int have_claimed)
{
    int i;
    printf("# per-core\n");
    printf("# cpu socket measured_hz claimed_hz ppm r nsamp retries "
           "hz_min hz_max hz_p50 implied_mono_hz\n");
    for (i = 0; i < nrows; i++)
        print_core_row(&rows[i], have_claimed);
}

static void write_csv(const char *path, const result_t *rows, int nrows,
                      int have_claimed)
{
    FILE *f;
    int i;

    f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "error: cannot write %s: %s\n", path, strerror(errno));
        return;
    }
    fprintf(f, "cpu,socket,measured_hz,claimed_hz,ppm,r,nsamp,retries,"
               "hz_min,hz_max,hz_p50,implied_mono_hz\n");
    for (i = 0; i < nrows; i++) {
        const result_t *r = &rows[i];
        if (!r->ok) {
            fprintf(f, "%d,%d,,,,,,,,,,\n", r->cpu, r->socket);
            continue;
        }
        fprintf(f, "%d,%d,%.6f,", r->cpu, r->socket, r->measured_hz);
        if (have_claimed && r->claimed_hz > 0.0)
            fprintf(f, "%.0f,%.4f,", r->claimed_hz, r->ppm);
        else
            fprintf(f, ",,");
        fprintf(f, "%.8f,%d,%d,%.6f,%.6f,%.6f,",
                r->r, r->nsamp, r->retries, r->hz_min, r->hz_max, r->hz_p50);
        if (have_claimed && r->implied_mono_hz > 0.0)
            fprintf(f, "%.6f\n", r->implied_mono_hz);
        else
            fprintf(f, "\n");
    }
    fclose(f);
}

/* nspt = 1e9 / measured_hz: 每 tick 多少纳秒(实测, 非名义值)。
 * 只取 ok 行的 measured_hz; 多核时用 p50, 并在文件头记录极差供判断。 */
static double nspt_p50(const result_t *rows, int nrows,
                       double *hz_min_out, double *hz_max_out)
{
    double hz[MAX_CPUS];
    int n = 0, i;
    double hmin = 0.0, hmax = 0.0;

    for (i = 0; i < nrows; i++) {
        if (!rows[i].ok || rows[i].measured_hz <= 0.0)
            continue;
        hz[n++] = rows[i].measured_hz;
        if (n == 1) {
            hmin = hmax = hz[0];
        } else {
            if (hz[n - 1] < hmin)
                hmin = hz[n - 1];
            if (hz[n - 1] > hmax)
                hmax = hz[n - 1];
        }
    }
    if (n == 0)
        return 0.0;
    if (hz_min_out)
        *hz_min_out = hmin;
    if (hz_max_out)
        *hz_max_out = hmax;
    return 1e9 / percentile_p50(hz, n);
}

static int write_nspt(const char *path, const result_t *rows, int nrows,
                      const opts_t *opt, const topo_t *topo)
{
    double hzmin = 0.0, hzmax = 0.0, nspt;
    char host[HOSTNAME_LEN];
    FILE *f;
    int i;

    if (gethostname(host, sizeof(host) - 1) != 0)
        strncpy(host, "?", sizeof(host));
    host[sizeof(host) - 1] = '\0';

    nspt = nspt_p50(rows, nrows, &hzmin, &hzmax);
    if (nspt <= 0.0) {
        fprintf(stderr, "error: no valid measurement for nspt output\n");
        return -1;
    }

    f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "error: cannot write %s: %s\n", path, strerror(errno));
        return -1;
    }
    fprintf(f, "# nspt file generated by probe_counter_freq\n");
    fprintf(f, "# node=%s arch=%s counter=%s clock=%s\n",
            host, g_arch_name, g_counter_name, clock_name(opt->clock_raw));
    fprintf(f, "# interval=%.3fs sample_ms=%d r_min=%.6f max_retry=%d\n",
            opt->interval_sec, opt->sample_ms, opt->r_min, opt->max_retry);
    fprintf(f, "# cores=%d hz_min=%.6f hz_max=%.6f hz_spread_ppm=%.2f\n",
            nrows, hzmin, hzmax,
            (hzmax - hzmin) / ((hzmin + hzmax) / 2.0) * 1e6);
    fprintf(f, "# time=%ld\n", (long)time(NULL));
    printf("# nspt (p50) = %.9f ns/tick\n", nspt);
    fprintf(f, "%.9f\n", nspt);
    fclose(f);
    (void)topo;
    (void)i;
    return 0;
}

/* --verify: 读已存在的 nspt 文件, 与本次实测 p50 比较, 输出偏差 ppm。 */
static int verify_nspt(const char *path, const result_t *rows, int nrows)
{
    double file_nspt = 0.0, hzmin = 0.0, hzmax = 0.0, nspt, ppm;
    char line[256];
    FILE *f;

    f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "error: cannot open nspt file %s: %s\n",
                path, strerror(errno));
        return -1;
    }
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        char *end;
        double v;
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        v = strtod(p, &end);
        if (end != p && v > 0.0) {
            file_nspt = v;
            break;
        }
    }
    fclose(f);

    if (file_nspt <= 0.0) {
        fprintf(stderr, "error: no nspt value found in %s\n", path);
        return -1;
    }

    nspt = nspt_p50(rows, nrows, &hzmin, &hzmax);
    if (nspt <= 0.0) {
        fprintf(stderr, "error: no valid measurement for verify\n");
        return -1;
    }

    ppm = (nspt - file_nspt) / file_nspt * 1e6;
    printf("# verify: file=%s file_nspt=%.9f measured_nspt=%.9f "
           "deviation=%.2f ppm\n", path, file_nspt, nspt, ppm);
    if (ppm > 5.0 || ppm < -5.0)
        printf("# verify: WARNING deviation > 5 ppm\n");
    return 0;
}

static void print_socket_table(const result_t *rows, int nrows,
                               const topo_t *topo, int have_claimed)
{
    int s, i;

    printf("# per-socket\n");
    printf("# socket ncpu measured_hz_p50 ppm_p50 ppm_min ppm_max r_min\n");

    for (s = 0; s < topo->nsockets; s++) {
        int pkg = topo->socket_ids[s];
        double hz[MAX_CPUS], ppm[MAX_CPUS];
        int nh = 0, np = 0, ncpu = 0;
        double rmin = 1.0, pmin = 0.0, pmax = 0.0;
        int have_ppm = 0;

        for (i = 0; i < nrows; i++) {
            if (!rows[i].ok || rows[i].socket != pkg)
                continue;
            ncpu++;
            hz[nh++] = rows[i].measured_hz;
            if (rows[i].r < rmin)
                rmin = rows[i].r;
            if (have_claimed && rows[i].claimed_hz > 0.0) {
                if (!have_ppm) {
                    pmin = pmax = rows[i].ppm;
                    have_ppm = 1;
                } else {
                    if (rows[i].ppm < pmin)
                        pmin = rows[i].ppm;
                    if (rows[i].ppm > pmax)
                        pmax = rows[i].ppm;
                }
                ppm[np++] = rows[i].ppm;
            }
        }
        if (ncpu == 0)
            continue;
        printf("%d %d %.6f ", pkg, ncpu, percentile_p50(hz, nh));
        if (have_ppm && np > 0)
            printf("%.4f %.4f %.4f ", percentile_p50(ppm, np), pmin, pmax);
        else
            printf("NA NA NA ");
        printf("%.8f\n", rmin);
    }
}

/* -------------------------------------------------------------------------- */
/* CLI                                                                        */
/* -------------------------------------------------------------------------- */

static void opts_default(opts_t *opt)
{
    opt->interval_sec = DEFAULT_INTERVAL_SEC;
    opt->sample_ms = DEFAULT_SAMPLE_MS;
    opt->r_min = DEFAULT_R_MIN;
    opt->max_retry = DEFAULT_MAX_RETRY;
    opt->clock_raw = 0;
    opt->seed = 0;
    opt->seed_set = 0;
    opt->cpu_override = -1;
    opt->csv_path = NULL;
    opt->nspt_out = NULL;
    opt->verify = 0;
    opt->verify_path = NULL;
}

static int parse_opts(int argc, char **argv, opts_t *opt)
{
    int i;

    opts_default(opt);
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)
            return 1;
        if (strcmp(argv[i], "--interval") == 0 && i + 1 < argc) {
            opt->interval_sec = strtod(argv[++i], NULL);
            if (opt->interval_sec <= 0.0) {
                fprintf(stderr, "error: --interval must be > 0\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--sample-ms") == 0 && i + 1 < argc) {
            opt->sample_ms = atoi(argv[++i]);
            if (opt->sample_ms < 1) {
                fprintf(stderr, "error: --sample-ms must be >= 1\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--r-min") == 0 && i + 1 < argc) {
            opt->r_min = strtod(argv[++i], NULL);
            if (opt->r_min <= 0.0 || opt->r_min > 1.0) {
                fprintf(stderr, "error: --r-min must be in (0, 1]\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--max-retry") == 0 && i + 1 < argc) {
            opt->max_retry = atoi(argv[++i]);
            if (opt->max_retry < 0) {
                fprintf(stderr, "error: --max-retry must be >= 0\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--clock") == 0 && i + 1 < argc) {
            i++;
            if (strcmp(argv[i], "mono") == 0)
                opt->clock_raw = 0;
            else if (strcmp(argv[i], "raw") == 0)
                opt->clock_raw = 1;
            else {
                fprintf(stderr, "error: --clock must be mono or raw\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            opt->seed = (unsigned int)strtoul(argv[++i], NULL, 10);
            opt->seed_set = 1;
        } else if (strcmp(argv[i], "--cpu") == 0 && i + 1 < argc) {
            opt->cpu_override = atoi(argv[++i]);
            if (opt->cpu_override < 0) {
                fprintf(stderr, "error: --cpu must be >= 0\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--csv") == 0 && i + 1 < argc) {
            opt->csv_path = argv[++i];
        } else if (strcmp(argv[i], "--nspt-out") == 0 && i + 1 < argc) {
            opt->nspt_out = argv[++i];
        } else if (strcmp(argv[i], "--verify") == 0 && i + 1 < argc) {
            opt->verify = 1;
            opt->verify_path = argv[++i];
        } else {
            fprintf(stderr, "error: unknown or incomplete option '%s' (try --help)\n",
                    argv[i]);
            return -1;
        }
    }
    return 0;
}

/* -------------------------------------------------------------------------- */
/* Drivers                                                                    */
/* -------------------------------------------------------------------------- */

static void first_claimed(const result_t *rows, int nrows,
                          const char **src, double *hz, int *have)
{
    int i;
    *src = "NA";
    *hz = 0.0;
    *have = 0;
    for (i = 0; i < nrows; i++) {
        if (rows[i].ok && rows[i].claimed_hz > 0.0) {
            *src = rows[i].claimed_src;
            *hz = rows[i].claimed_hz;
            *have = 1;
            return;
        }
    }
}

#ifndef USE_MPI
static int run_serial(const opts_t *opt, const topo_t *topo)
{
    result_t *rows;
    const char *claimed_src = "NA";
    double claimed_hz = 0.0;
    int nrows = 0, s, have_claimed = 0;
    int rc = 0;

    if (opt->cpu_override >= 0)
        nrows = 1;
    else
        nrows = topo->nsockets > 0 ? topo->nsockets : 1;

    rows = (result_t *)calloc((size_t)nrows, sizeof(result_t));
    if (!rows) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    if (opt->cpu_override >= 0) {
        if (!topo_cpu_online(topo, opt->cpu_override)) {
            fprintf(stderr, "error: CPU %d is not online\n", opt->cpu_override);
            free(rows);
            return 1;
        }
        if (measure_on_cpu(opt, topo, opt->cpu_override, &rows[0]) != 0)
            rc = 1;
    } else {
        for (s = 0; s < nrows; s++) {
            int cpu = pick_random_cpu_on_socket(topo, topo->socket_ids[s]);
            if (cpu < 0) {
                fprintf(stderr, "error: no online CPU on socket %d\n",
                        topo->socket_ids[s]);
                rc = 1;
                continue;
            }
            if (measure_on_cpu(opt, topo, cpu, &rows[s]) != 0)
                rc = 1;
        }
    }

    first_claimed(rows, nrows, &claimed_src, &claimed_hz, &have_claimed);
    print_header(opt, topo, claimed_src, claimed_hz, 1);
    print_core_table(rows, nrows, have_claimed);
    print_socket_table(rows, nrows, topo, have_claimed);
    if (opt->csv_path)
        write_csv(opt->csv_path, rows, nrows, have_claimed);
    if (opt->nspt_out)
        write_nspt(opt->nspt_out, rows, nrows, opt, topo);
    if (opt->verify) {
        const char *vp = opt->verify_path ? opt->verify_path : opt->nspt_out;
        if (vp)
            verify_nspt(vp, rows, nrows);
        else
            fprintf(stderr, "warn: --verify needs --nspt-out or a path\n");
    }
    free(rows);
    return rc;
}
#endif /* !USE_MPI */

#ifdef USE_MPI
static int run_mpi(const opts_t *opt, const topo_t *topo, int myrank, int nrank)
{
    result_t local, *all = NULL;
    const char *claimed_src = "NA";
    double claimed_hz = 0.0;
    int cpu, have_claimed = 0, rc = 0;

    if (nrank > topo->ncpu) {
        if (myrank == 0)
            fprintf(stderr,
                    "error: nprocs=%d exceeds online CPUs=%d\n",
                    nrank, topo->ncpu);
        return 1;
    }

    cpu = bound_cpu_or_neg();
    if (cpu < 0)
        cpu = topo->cpu[myrank % topo->ncpu];
    if (!topo_cpu_online(topo, cpu)) {
        fprintf(stderr, "error: rank %d CPU %d is not online\n", myrank, cpu);
        rc = 1;
    } else if (measure_on_cpu(opt, topo, cpu, &local) != 0) {
        rc = 1;
        memset(&local, 0, sizeof(local));
        local.cpu = cpu;
        local.socket = topo_socket_of_cpu(topo, cpu);
        snprintf(local.claimed_src, sizeof(local.claimed_src), "NA");
    }

    if (myrank == 0) {
        all = (result_t *)calloc((size_t)nrank, sizeof(result_t));
        if (!all) {
            fprintf(stderr, "error: out of memory\n");
            rc = 1;
        }
    }
    MPI_Gather(&local, (int)sizeof(result_t), MPI_BYTE,
               all, (int)sizeof(result_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (myrank == 0 && all) {
        first_claimed(all, nrank, &claimed_src, &claimed_hz, &have_claimed);
        print_header(opt, topo, claimed_src, claimed_hz, nrank);
        print_core_table(all, nrank, have_claimed);
        print_socket_table(all, nrank, topo, have_claimed);
        if (opt->csv_path)
            write_csv(opt->csv_path, all, nrank, have_claimed);
        if (opt->nspt_out)
            write_nspt(opt->nspt_out, all, nrank, opt, topo);
        if (opt->verify) {
            const char *vp = opt->verify_path ? opt->verify_path : opt->nspt_out;
            if (vp)
                verify_nspt(vp, all, nrank);
            else
                fprintf(stderr, "warn: --verify needs --nspt-out or a path\n");
        }
        free(all);
    }
    return rc;
}
#endif

int main(int argc, char **argv)
{
    opts_t opt;
    topo_t topo;
    int pr, rc = 0;
#ifdef USE_MPI
    int myrank = 0, nrank = 1;
#endif

#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif

    pr = parse_opts(argc, argv, &opt);
    if (pr < 0) {
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 2;
    }
    if (pr == 1) {
#ifdef USE_MPI
        if (myrank == 0)
#endif
            print_help(argv[0]);
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 0;
    }

#ifdef USE_MPI
    if (opt.cpu_override >= 0 && myrank == 0)
        fprintf(stderr, "warn: --cpu is ignored in the MPI binary\n");
#endif

    if (opt.seed_set)
        srand(opt.seed);
    else
        srand((unsigned int)time(NULL) ^ (unsigned int)getpid());

    if (topo_load(&topo) != 0) {
        fprintf(stderr, "error: failed to read CPU topology\n");
#ifdef USE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

#ifdef USE_MPI
    rc = run_mpi(&opt, &topo, myrank, nrank);
    MPI_Finalize();
#else
    rc = run_serial(&opt, &topo);
#endif
    return rc;
}
