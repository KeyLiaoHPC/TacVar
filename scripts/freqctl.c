/*
 * freqctl.c — adjtimex 时钟频率控制小工具（show / zero / nano）
 *
 * 用途：关闭 NTP 服务后，把内核时钟的驯服残留（freq/offset/status）清零，
 *       避免 CLOCK_MONOTONIC / CLOCK_REALTIME 相对 tick 时钟继续漂移。
 * 用法：
 *   freqctl          # 显示当前 adjtimex 状态
 *   freqctl zero     # 清零 freq + offset + status（保留 STA_NANO 位）（需要 root/sudo）
 *   freqctl pll-nano # 置 status = STA_PLL|STA_NANO（ARM 内核单独设 NANO 会被忽略）
 *   freqctl freq     # 只清零 freq
 * 编译：cc -O2 -std=c99 freqctl.c -o freqctl
 *
 * 已知坑（c920bn1, 2026-08-16 实测）：adjtimex status 清 0 后，ARM 内核走
 * USEC tick 路径，CLOCK_MONOTONIC/REALTIME 相对真实时间快约 100 ppm
 * （cntvct vs MONOTONIC_RAW = -0.43 ppm，vs MONOTONIC = -100.4 ppm）。
 * 保持 STA_NANO 位（0x2000）可回到正常速率。
 */
#include <stdio.h>
#include <string.h>
#include <sys/timex.h>

static int show(struct timex *tx)
{
    tx->modes = 0;   /* 只读查询 */
    if (adjtimex(tx) < 0) {
        perror("adjtimex");
        return 1;
    }
    printf("freq=%ld offset=%ld status=0x%x maxerror=%ld esterror=%ld\n",
           tx->freq, tx->offset, tx->status, tx->maxerror, tx->esterror);
    return 0;
}

int main(int argc, char **argv)
{
    struct timex tx;
    memset(&tx, 0, sizeof(tx));
    if (argc > 1 && strcmp(argv[1], "zero") == 0) {
        /* 保留 STA_NANO：清掉会触发 usec tick 路径的 -100ppm 伪影（ARM 实测） */
        if (adjtimex(&tx) < 0) { perror("adjtimex"); return 1; }
        tx.modes = ADJ_FREQUENCY | ADJ_OFFSET | ADJ_STATUS;
        tx.freq = 0;
        tx.offset = 0;
        tx.status = tx.status & STA_NANO;
    } else if (argc > 1 && strcmp(argv[1], "nano") == 0) {
        if (adjtimex(&tx) < 0) { perror("adjtimex"); return 1; }
        tx.modes = ADJ_STATUS;
        tx.status = STA_NANO;
    } else if (argc > 1 && strcmp(argv[1], "pll-nano") == 0) {
        if (adjtimex(&tx) < 0) { perror("adjtimex"); return 1; }
        tx.modes = ADJ_STATUS;
        tx.status = STA_PLL | STA_NANO;
    } else if (argc > 1 && strcmp(argv[1], "freq") == 0) {
        tx.modes = ADJ_FREQUENCY;
        tx.freq = 0;
    }
    return show(&tx);
}