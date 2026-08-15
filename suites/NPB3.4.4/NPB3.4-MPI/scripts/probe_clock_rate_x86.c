/*
 * probe_clock_rate_x86.c — x86 port of probe_clock_rate.c (c920bn1 aarch64)
 * 目的: 在 x86-64 节点上量化 TSC(rdtscp) 与 clock_gettime / papi_get_real_nsec
 *       之间的速率差(ppm)与漂移稳定性, 复刻 aarch64 cntvct 79ppm 判定方法。
 * 与 aarch64 版的差异:
 *   - 反汇编计数器: rdtscp (无 lfence, 与 gauge_tf_insitu.h 热路径一致), 记录 AUX
 *   - tick->ns 用 CPUID 0x15/0x16 名义频率 (印刷 TSC MHz), 不再固定 NSTP=10
 *   - 时钟集: CLOCK_MONOTONIC / CLOCK_MONOTONIC_RAW / CLOCK_REALTIME /
 *     PAPI_get_real_nsec (dlopen /home/hpckey/01-App/papi/lib/libpapi.so.7.2,
 *     失败则 papi 列记 -1)
 *
 * 用法:
 *   probe_clock_rate_x86 --abs [N]
 *       打印 N 轮 (默认 10) 括号读: 每轮 t0 m t1 r t2 k t3 p t4
 *       输出: aux0 t0 t1 t2 t3 t4 mono raw realtime papi (整数 ns/ticks)
 *   probe_clock_rate_x86 --gauge NG NITER [CLOCK]
 *       subtraction gauge, CLOCK=mono|raw|realtime|papi (默认 mono)
 *       ratio = tick_ns / clock_ns, ppm = (1-mean)*1e6 (>0 => tick 慢)
 *   probe_clock_rate_x86 --gauge-sym NG NITER [CLOCK]   对称括号取均值
 *   probe_clock_rate_x86 --series NG NITER ROUNDS GAP [CLOCK]   分轮稳定性
 *   probe_clock_rate_x86 --long SEC
 *       每 10 s 一轮绝对括号读, 持续 SEC 秒, 供 Python 拟合漂移斜率
 *   probe_clock_rate_x86 --sweep CPU
 *       单核短测: 5 轮 abs + 1 次 gauge(513040,100), 打印一行 per-core 结果
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <math.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* ---------------- TSC 读取 ---------------- */
static inline uint64_t rd_tsc(unsigned *aux)
{
    unsigned hi, lo;
    __asm__ volatile("rdtscp" : "=a"(lo), "=d"(hi), "=c"(*aux) :: "memory");
    return ((uint64_t)hi << 32) | (uint64_t)lo;
}

/* ---------------- 时钟读取 ---------------- */
static inline uint64_t rd_clock(int which)
{
    struct timespec ts;
    switch (which) {
    case 0: clock_gettime(CLOCK_MONOTONIC, &ts); break;
    case 1: clock_gettime(CLOCK_MONOTONIC_RAW, &ts); break;
    default: clock_gettime(CLOCK_REALTIME, &ts); break;
    }
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/* PAPI dlopen: 静态 os vector 已含 get_real_nsec, 无需 library_init (源码核实) */
static void *papi_handle = NULL;
static long long (*papi_real_nsec)(void) = NULL;

static void papi_init(void)
{
    const char *paths[] = {
        getenv("TACVAR_PAPI_LIB"),
        "/home/hpckey/01-App/papi/lib/libpapi.so.7.2",
        "/home/hpckey/01-App/papi/lib/libpapi.so.7",
        "/home/hpckey/01-App/papi/lib/libpapi.so",
        NULL
    };
    for (int i = 0; i < (int)(sizeof(paths) / sizeof(paths[0])); i++) {
        if (!paths[i]) continue;
        papi_handle = dlopen(paths[i], RTLD_NOW | RTLD_LOCAL);
        if (papi_handle) break;
        fprintf(stderr, "  dlopen(%s) failed: %s\n", paths[i], dlerror());
    }
    if (papi_handle) {
        papi_real_nsec = (long long (*)(void))dlsym(papi_handle, "PAPI_get_real_nsec");
        if (!papi_real_nsec) {
            fprintf(stderr, "warn: PAPI_get_real_nsec not found in libpapi\n");
            dlclose(papi_handle);
            papi_handle = NULL;
        }
    } else {
        fprintf(stderr, "warn: libpapi not loadable, papi column = -1 (%s)\n", dlerror());
    }
}

static inline long long rd_papi(void)
{
    return papi_real_nsec ? papi_real_nsec() : -1LL;
}

/* ---------------- TSC 频率标定 (名义值) ---------------- */
static double g_ns_per_tick = 0.0;   /* 名义 ns/tick */
static double g_tsc_ghz = 0.0;

static void tsc_freq_print(void)
{
    unsigned eax, ebx, ecx, edx;
    __asm__ volatile("cpuid" : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) : "a"(0x15));
    if (ebx != 0) {
        unsigned long long tsc_hz = (unsigned long long)ecx * ebx / eax;
        g_ns_per_tick = 1e9 / (double)tsc_hz;
        g_tsc_ghz = (double)tsc_hz / 1e9;
        printf("TSC nominal: CPUID.15H = %llu Hz (%.6f GHz), ns/tick=%.9f\n",
               (unsigned long long)tsc_hz, g_tsc_ghz, g_ns_per_tick);
        return;
    }
    __asm__ volatile("cpuid" : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) : "a"(0x16));
    if (eax != 0) {
        double mhz = (double)eax;
        g_ns_per_tick = 1e3 / mhz;
        g_tsc_ghz = mhz / 1e3;
        printf("TSC nominal: CPUID.16H = %.0f MHz (%.6f GHz), ns/tick=%.9f\n",
               mhz, g_tsc_ghz, g_ns_per_tick);
        return;
    }
    /* 最后兜底: 500 ms 窗口对 CLOCK_MONOTONIC 实测 */
    unsigned aux;
    uint64_t t0 = rd_tsc(&aux);
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    uint64_t g0 = (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
    usleep(500000);
    clock_gettime(CLOCK_MONOTONIC, &ts);
    uint64_t g1 = (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
    uint64_t t1 = rd_tsc(&aux);
    double hz = (double)(t1 - t0) / ((double)(g1 - g0) * 1e-9);
    g_ns_per_tick = 1e9 / hz;
    g_tsc_ghz = hz / 1e9;
    printf("TSC nominal: fallback measured = %.6f GHz, ns/tick=%.9f\n",
           g_tsc_ghz, g_ns_per_tick);
}

static inline double tick_ns_delta(uint64_t b, uint64_t e)
{
    return (double)(e - b) * g_ns_per_tick;
}

/* ---------------- 工具 ---------------- */
static void pin_cpu(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (sched_setaffinity(0, sizeof(set), &set) != 0)
        fprintf(stderr, "warn: sched_setaffinity(%d) failed: %s\n", cpu,
                strerror(errno));
}

static inline void gauge_loop(uint64_t nsteps, uint64_t seed)
{
    register uint64_t ra = seed;
    for (uint64_t i = 0; i < nsteps; i++) {
        /* 与真实插桩同构: 逐条 sub, 防 -O2 折叠 */
        __asm__ volatile("subq $1, %0" : "+r"(ra));
    }
    (void)ra;
}

static const char *clock_name(int which)
{
    switch (which) {
    case 0: return "monotonic";
    case 1: return "monotonic_raw";
    case 2: return "realtime";
    case 3: return "papi";
    default: return "?";
    }
}

/* ---------------- 模式实现 ---------------- */

/* --abs [N]: 括号读, 每轮 t0 m t1 r t2 k t3 p t4 */
static void run_abs(int n)
{
    printf("# abs rounds=%d  (t0..t4=rdtscp ticks, m=mono, r=raw, k=realtime, p=papi; ns)\n", n);
    for (int i = 0; i < n; i++) {
        unsigned a0, a1, a2, a3, a4;
        uint64_t t0 = rd_tsc(&a0);
        uint64_t m  = rd_clock(0);
        uint64_t t1 = rd_tsc(&a1);
        uint64_t r  = rd_clock(1);
        uint64_t t2 = rd_tsc(&a2);
        uint64_t k  = rd_clock(2);
        uint64_t t3 = rd_tsc(&a3);
        long long p = rd_papi();
        uint64_t t4 = rd_tsc(&a4);
        printf("%u %llu %llu %llu %llu %llu %llu %llu %llu %lld\n",
               a0,
               (unsigned long long)t0, (unsigned long long)t1,
               (unsigned long long)t2, (unsigned long long)t3,
               (unsigned long long)t4,
               (unsigned long long)m, (unsigned long long)r,
               (unsigned long long)k, (long long)p);
    }
}

/* --gauge NG NITER [CLOCK] */
static void run_gauge(uint64_t ng, uint64_t niters, int which, int quiet)
{
    uint64_t seed = 0x123456789abcdef0ULL;
    double sum_r = 0.0, sum_r2 = 0.0, rmin = 1e18, rmax = 0.0;
    uint64_t n = 0;
    for (uint64_t i = 0; i < niters; i++) {
        gauge_loop(64, seed);
        unsigned a0, a1;
        uint64_t c0 = rd_tsc(&a0);
        uint64_t g0 = rd_clock(which);
        if (which == 3) g0 = (uint64_t)rd_papi();
        gauge_loop(ng, seed);
        uint64_t g1 = rd_clock(which);
        if (which == 3) g1 = (uint64_t)rd_papi();
        uint64_t c1 = rd_tsc(&a1);
        int64_t dg = (int64_t)(g1 - g0);
        if (dg <= 0) continue;
        double r = tick_ns_delta(c0, c1) / (double)dg;   /* tick-ns per clock-ns */
        sum_r += r;
        sum_r2 += r * r;
        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
        n++;
        if (!quiet && i < 5)
            printf("  gauge[%llu] dg=%lld ns  ratio=%.6f\n",
                   (unsigned long long)i, (long long)dg, r);
    }
    if (n == 0) { printf("  no valid gauges\n"); return; }
    double mean = sum_r / n;
    double var = sum_r2 / n - mean * mean;
    double sd = (var > 0) ? sqrt(var) : 0.0;
    printf("  clock=%-8s ng=%-10llu n=%llu  ratio mean=%.6f  sd=%.6f  min=%.6f  max=%.6f  (1-mean)=%.1f ppm\n",
           clock_name(which), (unsigned long long)ng, (unsigned long long)n,
           mean, sd, rmin, rmax, (1.0 - mean) * 1e6);
}

/* --gauge-sym NG NITER [CLOCK] */
static void run_gauge_sym(uint64_t ng, uint64_t niters, int which)
{
    uint64_t seed = 0x123456789abcdef0ULL;
    double sumA = 0.0, sumB = 0.0;
    uint64_t nA = 0, nB = 0;
    for (uint64_t i = 0; i < niters; i++) {
        gauge_loop(64, seed);
        unsigned a0, a1;
        uint64_t g0 = rd_clock(which), c0;
        if (which == 3) g0 = (uint64_t)rd_papi();
        c0 = rd_tsc(&a0);
        gauge_loop(ng, seed);
        uint64_t c1 = rd_tsc(&a1), g1;
        if (which == 3) g1 = (uint64_t)rd_papi();
        else g1 = rd_clock(which);
        int64_t dg = (int64_t)(g1 - g0);
        if (dg > 0) { sumB += tick_ns_delta(c0, c1) / (double)dg; nB++; }
        /* 模式 A */
        uint64_t t0 = rd_tsc(&a0);
        uint64_t h0 = rd_clock(which);
        if (which == 3) h0 = (uint64_t)rd_papi();
        gauge_loop(ng, seed);
        uint64_t h1 = rd_clock(which);
        if (which == 3) h1 = (uint64_t)rd_papi();
        uint64_t t1 = rd_tsc(&a1);
        int64_t dh = (int64_t)(h1 - h0);
        if (dh > 0) { sumA += tick_ns_delta(t0, t1) / (double)dh; nA++; }
    }
    double mA = sumA / nA, mB = sumB / nB;
    printf("  clock=%-8s sym-ng=%-10llu A=%.6f B=%.6f  avg=%.6f  -> 偏差 %.1f ppm\n",
           clock_name(which), (unsigned long long)ng, mA, mB, (mA + mB) / 2.0,
           (1.0 - (mA + mB) / 2.0) * 1e6);
}

/* --series NG NITER ROUNDS GAP [CLOCK] */
static void run_series(uint64_t ng, uint64_t niters, int rounds, int gap, int which)
{
    printf("series: ng=%llu niters=%llu rounds=%d gap=%ds clock=%s\n",
           (unsigned long long)ng, (unsigned long long)niters, rounds, gap,
           clock_name(which));
    for (int r = 0; r < rounds; r++) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        printf("round %d @ wall %lld.%09ld s:\n", r,
               (long long)ts.tv_sec, ts.tv_nsec);
        run_gauge(ng, niters, which, 1);
        if (r + 1 < rounds)
            sleep(gap);
    }
}

/* --long SEC: 每 10 s 一轮 --abs 格式, 供拟合 */
static void run_long(int seconds)
{
    printf("# long window %d s, one bracket round per 10 s: wall_ns aux t0 t1 t2 t3 t4 m r k p\n",
           seconds);
    unsigned a0, a1, a2, a3, a4;
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    uint64_t wall = (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
    int rounds = seconds / 10;
    for (int i = 0; i < rounds; i++) {
        clock_gettime(CLOCK_MONOTONIC, &ts);
        wall = (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
        uint64_t t0 = rd_tsc(&a0);
        uint64_t m  = rd_clock(0);
        uint64_t t1 = rd_tsc(&a1);
        uint64_t r  = rd_clock(1);
        uint64_t t2 = rd_tsc(&a2);
        uint64_t k  = rd_clock(2);
        uint64_t t3 = rd_tsc(&a3);
        long long p = rd_papi();
        uint64_t t4 = rd_tsc(&a4);
        printf("%llu %u %llu %llu %llu %llu %llu %llu %llu %llu %lld\n",
               (unsigned long long)wall, a0,
               (unsigned long long)t0, (unsigned long long)t1,
               (unsigned long long)t2, (unsigned long long)t3,
               (unsigned long long)t4,
               (unsigned long long)m, (unsigned long long)r,
               (unsigned long long)k, (long long)p);
        fflush(stdout);
        if (i + 1 < rounds)
            sleep(10);
    }
}

/* --sweep CPU: 单核短测 */
static void run_sweep(int cpu)
{
    pin_cpu(cpu);
    gauge_loop(1000, 0x1234);         /* warmup */
    run_abs(5);
    run_gauge(513040, 100, 0, 1);      /* mono */
    run_gauge(513040, 100, 1, 1);      /* raw */
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr,
                "usage: %s --abs [N] | --gauge NG NITER [CLOCK] | "
                "--gauge-sym NG NITER [CLOCK] | --series NG NITER ROUNDS GAP [CLOCK] | "
                "--long SEC | --sweep CPU\n",
                argv[0]);
        return 2;
    }
    papi_init();
    tsc_freq_print();
    char hname[256] = "?";
    if (gethostname(hname, sizeof(hname) - 1) != 0)
        strcpy(hname, "?");
    printf("# node: %s  pid=%d\n", hname, getpid());

    if (strcmp(argv[1], "--abs") == 0) {
        int n = (argc >= 3) ? atoi(argv[2]) : 10;
        if (n <= 0) n = 10;
        pin_cpu(0);
        run_abs(n);
        return 0;
    }
    if (strcmp(argv[1], "--gauge") == 0 && argc >= 4) {
        uint64_t ng = strtoull(argv[2], NULL, 0);
        uint64_t ni = strtoull(argv[3], NULL, 0);
        int which = (argc >= 5) ? atoi(argv[4]) : 0;
        pin_cpu(0);
        printf("gauge-based, ng=%llu niters=%llu clock=%s:\n",
               (unsigned long long)ng, (unsigned long long)ni, clock_name(which));
        run_gauge(ng, ni, which, 0);
        return 0;
    }
    if (strcmp(argv[1], "--gauge-sym") == 0 && argc >= 4) {
        uint64_t ng = strtoull(argv[2], NULL, 0);
        uint64_t ni = strtoull(argv[3], NULL, 0);
        int which = (argc >= 5) ? atoi(argv[4]) : 0;
        pin_cpu(0);
        printf("symmetric-bracket gauge, ng=%llu niters=%llu clock=%s:\n",
               (unsigned long long)ng, (unsigned long long)ni, clock_name(which));
        run_gauge_sym(ng, ni, which);
        return 0;
    }
    if (strcmp(argv[1], "--series") == 0 && argc >= 6) {
        uint64_t ng = strtoull(argv[2], NULL, 0);
        uint64_t ni = strtoull(argv[3], NULL, 0);
        int rounds = atoi(argv[4]);
        int gap = atoi(argv[5]);
        int which = (argc >= 7) ? atoi(argv[6]) : 0;
        pin_cpu(0);
        run_series(ng, ni, rounds, gap, which);
        return 0;
    }
    if (strcmp(argv[1], "--long") == 0 && argc >= 3) {
        int sec = atoi(argv[2]);
        if (sec < 20) sec = 120;
        pin_cpu(0);
        run_long(sec);
        return 0;
    }
    if (strcmp(argv[1], "--sweep") == 0 && argc >= 3) {
        int cpu = atoi(argv[2]);
        printf("===== sweep cpu=%d =====\n", cpu);
        fflush(stdout);
        run_sweep(cpu);
        return 0;
    }
    fprintf(stderr, "bad args\n");
    return 2;
}