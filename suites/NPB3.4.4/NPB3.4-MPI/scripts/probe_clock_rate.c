/*
 * probe_clock_rate.c — cntvct_el0 vs clock_gettime(CLOCK_MONOTONIC) 速率/漂移探针
 * 用于验证 c920bn1 上两时钟 79ppm 偏差是否稳定（速率恒定 vs 时变漂移）。
 *
 * 用法:
 *   probe_clock_rate --abs              打印 10 组 (clock_ns, cntvct_ticks) 绝对读数
 *   probe_clock_rate --gauge NG NITER   跑 NITER 次 ng 步 subtraction gauge，
 *                                       每轮同时记录两时钟，输出 per-gauge 比率统计
 *   probe_clock_rate --series NG NITER ROUNDS GAP_S
 *                                       分 ROUNDS 轮测量，每轮后 sleep GAP_S，
 *                                       观察比率随墙钟时间的稳定性
 *
 * 说明:
 *   - gauge 循环与真实插桩同为"减法计步"结构（rdtscp/cntvct 无屏障读）。
 *   - 每轮 gauge 的读序: c0(cntvct) g0(clock) [gauge] g1(clock) c1(cntvct)
 *     比率 = (c1-c0)*NSTP / (g1-g0)，NSTP 对 cntvct 固定 10。
 *   - 绑定单核 (sched_setaffinity) 避免迁移噪声。
 */
#define _GNU_SOURCE
#include <errno.h>
#include <math.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define NSTP 10

static inline uint64_t rd_cntvct(void)
{
    uint64_t v;
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(v));
    return v;
}

static inline uint64_t rd_clock(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static inline void gauge_loop(uint64_t nsteps, uint64_t seed)
{
    register uint64_t ra = seed;
    for (uint64_t i = 0; i < nsteps; i++) {
        /* 与真实插桩同构：逐条 sub 指令，防 -O2 折叠/删除循环 */
        __asm__ volatile("sub %0, %0, #1" : "+r"(ra));
    }
    (void)ra;
}

static void pin_cpu(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (sched_setaffinity(0, sizeof(set), &set) != 0)
        fprintf(stderr, "warn: sched_setaffinity(%d) failed: %s\n", cpu,
                strerror(errno));
}

static void run_abs(void)
{
    for (int i = 0; i < 10; i++) {
        uint64_t c = rd_cntvct();
        uint64_t g = rd_clock();
        /* 再读一次 cntvct，取两次的平均刻数，降低读取顺序噪声 */
        uint64_t c2 = rd_cntvct();
        printf("%llu %llu\n", (unsigned long long)g,
               (unsigned long long)((c + c2) / 2));
    }
}

static void run_gauge(uint64_t ng, uint64_t niters, int quiet)
{
    uint64_t seed = 0x123456789abcdef0ULL;
    double sum_r = 0.0, sum_r2 = 0.0, rmin = 1e18, rmax = 0.0;
    uint64_t n = 0;
    for (uint64_t i = 0; i < niters; i++) {
        gauge_loop(64, seed);            /* warmup */
        uint64_t c0 = rd_cntvct();
        uint64_t g0 = rd_clock();
        gauge_loop(ng, seed);
        uint64_t g1 = rd_clock();
        uint64_t c1 = rd_cntvct();
        int64_t dg = (int64_t)(g1 - g0);
        int64_t dc = (int64_t)(c1 - c0);
        if (dg <= 0 || dc <= 0)
            continue;
        double r = ((double)dc * NSTP) / (double)dg;   /* cntvct-ns per clock-ns */
        sum_r += r;
        sum_r2 += r * r;
        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
        n++;
        if (!quiet && i < 5)
            printf("  gauge[%llu] dg=%lld ns  dc=%lld ticks  ratio=%.6f\n",
                   (unsigned long long)i, (long long)dg, (long long)dc, r);
    }
    if (n == 0) {
        printf("  no valid gauges\n");
        return;
    }
    double mean = sum_r / n;
    double var = sum_r2 / n - mean * mean;
    double sd = (var > 0) ? sqrt(var) : 0.0;
    printf("  ng=%-10llu n=%llu  ratio mean=%.6f  sd=%.6f  min=%.6f  max=%.6f  (1-mean)=%.1f ppm\n",
           (unsigned long long)ng, (unsigned long long)n, mean, sd, rmin, rmax,
           (1.0 - mean) * 1e6);
}

/* 对称括号 gauge：交替两种读序并取均值，抵消单侧读钟间隙偏差 */
static void run_gauge_sym(uint64_t ng, uint64_t niters)
{
    uint64_t seed = 0x123456789abcdef0ULL;
    double sumA = 0.0, sumB = 0.0;
    uint64_t nA = 0, nB = 0;
    for (uint64_t i = 0; i < niters; i++) {
        gauge_loop(64, seed);
        if (i & 1) {                       /* 模式 B: g0 c0 |g| c1 g1 */
            uint64_t g0 = rd_clock(), c0 = rd_cntvct();
            gauge_loop(ng, seed);
            uint64_t c1 = rd_cntvct(), g1 = rd_clock();
            int64_t dg = (int64_t)(g1 - g0), dc = (int64_t)(c1 - c0);
            if (dg > 0 && dc > 0) { sumB += ((double)dc * NSTP) / dg; nB++; }
        } else {                           /* 模式 A: c0 g0 |g| g1 c1 */
            uint64_t c0 = rd_cntvct(), g0 = rd_clock();
            gauge_loop(ng, seed);
            uint64_t g1 = rd_clock(), c1 = rd_cntvct();
            int64_t dg = (int64_t)(g1 - g0), dc = (int64_t)(c1 - c0);
            if (dg > 0 && dc > 0) { sumA += ((double)dc * NSTP) / dg; nA++; }
        }
    }
    double mA = sumA / nA, mB = sumB / nB;
    printf("  sym-ng=%-10llu A=%.6f B=%.6f  avg=%.6f  -> 偏差 %.1f ppm\n",
           (unsigned long long)ng, mA, mB, (mA + mB) / 2.0,
           (1.0 - (mA + mB) / 2.0) * 1e6);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s --abs | --gauge NG NITER | --series NG NITER ROUNDS GAP_S\n", argv[0]);
        return 2;
    }
    pin_cpu(0);

    if (strcmp(argv[1], "--abs") == 0) {
        run_abs();
        return 0;
    }
    if (strcmp(argv[1], "--gauge") == 0 && argc >= 4) {
        uint64_t ng = strtoull(argv[2], NULL, 0);
        uint64_t ni = strtoull(argv[3], NULL, 0);
        printf("gauge-based, ng=%llu niters=%llu:\n",
               (unsigned long long)ng, (unsigned long long)ni);
        run_gauge(ng, ni, 0);
        return 0;
    }
    if (strcmp(argv[1], "--gauge-sym") == 0 && argc >= 4) {
        uint64_t ng = strtoull(argv[2], NULL, 0);
        uint64_t ni = strtoull(argv[3], NULL, 0);
        printf("symmetric-bracket gauge, ng=%llu niters=%llu:\n",
               (unsigned long long)ng, (unsigned long long)ni);
        run_gauge_sym(ng, ni);
        return 0;
    }
    if (strcmp(argv[1], "--series") == 0 && argc >= 6) {
        uint64_t ng = strtoull(argv[2], NULL, 0);
        uint64_t ni = strtoull(argv[3], NULL, 0);
        int rounds = atoi(argv[4]);
        int gap = atoi(argv[5]);
        printf("series: ng=%llu niters=%llu rounds=%d gap=%ds\n",
               (unsigned long long)ng, (unsigned long long)ni, rounds, gap);
        for (int r = 0; r < rounds; r++) {
            struct timespec ts;
            clock_gettime(CLOCK_MONOTONIC, &ts);
            printf("round %d @ wall %lld.%09ld s:\n", r,
                   (long long)ts.tv_sec, ts.tv_nsec);
            run_gauge(ng, ni, 1);
            if (r + 1 < rounds)
                sleep(gap);
        }
        return 0;
    }
    fprintf(stderr, "bad args\n");
    return 2;
}