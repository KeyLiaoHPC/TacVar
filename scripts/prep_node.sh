#!/usr/bin/env bash
# =============================================================================
# prep_node.sh — 测量节点环境准备 / 状态查看 / 恢复
#   1) 锁频：cpufreq 治理器改 performance，并把每个核的 scaling_min/max 钳到
#      LOCK_KHZ（默认 = scaling_max_freq；建议按节点勘察结果显式传入基频档）
#   2) NTP：停 chronyd / systemd-timesyncd，adjtimex 清零 freq+offset+status
#      （仅停服务会残留内核驯服频率，导致 CLOCK_MONOTONIC/REALTIME 相对 tick 漂移）
#      注意（c920bn1, 2026-08-16 实测）：openEuler 5.10 ARM 上 adjtimex 清零后
#      MONOTONIC/REALTIME 会出现 +100ppm 的 USEC-tick 伪影（cntvct vs RAW 正常）。
#      该节点改用"冻结"协议：先起 chronyd 驯服稳定，再 stop 服务且不做任何
#      adjtimex（冻结驯服态），再按冻结态实测 nspt。详见 freqctl.c 注释。
#   3) ssh：写入 /etc/ssh/sshd_config.d/99-allow-hpckey.conf（AllowUsers hpckey）
#      禁止其它账号新登录（不杀已有会话），收尾时 restore 删除并恢复
# 用法（在目标节点上以 hpckey 运行，sudo 免密）：
#   bash scripts/prep_node.sh setup   <NODE> [LOCK_KHZ]
#   bash scripts/prep_node.sh status  <NODE>
#   bash scripts/prep_node.sh restore <NODE>
# 原值记录在 /home/hpckey/.tacvar_prep/<NODE>/orig.env（restore 依此恢复）。
# =============================================================================
set -uo pipefail

ACTION="${1:-status}"
NODE="${2:?用法: prep_node.sh setup|status|restore NODE [LOCK_KHZ]}"
LOCK_KHZ="${3:-}"
REC_DIR="/home/hpckey/.tacvar_prep/${NODE}"
sudo -n true 2>/dev/null || { echo "sudo -n 不可用，退出"; exit 1; }

log() { echo "[$(date '+%F %T')] $*"; }

record_orig() {
  mkdir -p "$REC_DIR"
  SYS=/sys/devices/system/cpu/cpu0/cpufreq
  {
    echo "governor=$(cat "$SYS/scaling_governor" 2>/dev/null)"
    echo "scaling_min_freq=$(cat "$SYS/scaling_min_freq" 2>/dev/null)"
    echo "scaling_max_freq=$(cat "$SYS/scaling_max_freq" 2>/dev/null)"
    echo "scaling_driver=$(cat "$SYS/scaling_driver" 2>/dev/null)"
    echo "chronyd_active=$(systemctl is-active chronyd 2>/dev/null || echo inactive)"
    echo "timesyncd_active=$(systemctl is-active systemd-timesyncd 2>/dev/null || echo inactive)"
    echo "sshd_dropin_present=$([[ -f /etc/ssh/sshd_config.d/99-allow-hpckey.conf ]] && echo yes || echo no)"
  } > "$REC_DIR/orig.env"
  cpupower frequency-info > "$REC_DIR/freqinfo_orig.txt" 2>&1 || true
  chmod 600 "$REC_DIR/orig.env"
  log "原值已记录: $REC_DIR/orig.env"
}

sysfs_get() {  # sysfs_get cpuN field
  cat "/sys/devices/system/cpu/cpu$1/cpufreq/$2" 2>/dev/null
}

setup() {
  record_orig
  # 1) 锁频
  log "治理器 -> performance（全部核心）"
  sudo cpupower -c all frequency-set -g performance || log "WARN: cpupower -g performance 失败"
  if [[ -z "$LOCK_KHZ" ]]; then
    LOCK_KHZ=$(sysfs_get 0 scaling_max_freq)
    log "LOCK_KHZ 未指定，取 scaling_max_freq = $LOCK_KHZ kHz"
  fi
  ncpu=0
  for f in /sys/devices/system/cpu/cpu[0-9]*/cpufreq/scaling_max_freq; do
    echo "$LOCK_KHZ" | sudo tee "$f" > /dev/null || break
    ncpu=$((ncpu+1))
  done
  nmin=0
  for f in /sys/devices/system/cpu/cpu[0-9]*/cpufreq/scaling_min_freq; do
    echo "$LOCK_KHZ" | sudo tee "$f" > /dev/null || break
    nmin=$((nmin+1))
  done
  log "锁频完成: 全部 $ncpu 核 max/min = $LOCK_KHZ kHz"
  # 2) NTP 关闭 + 内核时钟驯服清零
  if systemctl is-active --quiet chronyd; then
    sudo systemctl stop chronyd && log "chronyd 已停止"
  fi
  if systemctl is-active --quiet systemd-timesyncd; then
    sudo systemctl stop systemd-timesyncd && log "systemd-timesyncd 已停止"
  fi
  FREQCTL=/home/hpckey/bin/freqctl
  if [[ ! -x "$FREQCTL" ]]; then
    mkdir -p /home/hpckey/bin
    cc -O2 -std=c99 "$(dirname "${BASH_SOURCE[0]}")/freqctl.c" -o "$FREQCTL" || { echo "freqctl 编译失败"; exit 1; }
    log "freqctl 编译完成"
  fi
  sudo "$FREQCTL" zero
  log "adjtimex 已清零（freq/offset/status），须静置约 150 秒让内核倍率收敛"
  # 3) ssh 禁新登录（仅 hpckey）
  echo 'AllowUsers hpckey' | sudo tee /etc/ssh/sshd_config.d/99-allow-hpckey.conf > /dev/null
  sudo systemctl reload ssh 2>/dev/null || sudo systemctl reload sshd 2>/dev/null || sudo kill -HUP "$(cat /run/sshd.pid 2>/dev/null)" 2>/dev/null
  sleep 1
  log "sshd AllowUsers 生效情况: $(sudo sshd -T 2>/dev/null | grep -i '^allowusers' || echo '无法确认')"
  log "setup 完成。请静置 ≥150 秒后再开始测量。"
}

status() {
  SYS=/sys/devices/system/cpu/cpu0/cpufreq
  echo "== 节点 $NODE 状态 =="
  echo "governor:   $(sysfs_get 0 scaling_governor)"
  echo "min/max:    $(sysfs_get 0 scaling_min_freq) / $(sysfs_get 0 scaling_max_freq) kHz"
  echo "cur:        $(sysfs_get 0 scaling_cur_freq) kHz"
  echo "load:       $(cat /proc/loadavg)"
  echo "chronyd:    $(systemctl is-active chronyd 2>/dev/null || echo inactive)"
  echo "timesyncd:  $(systemctl is-active systemd-timesyncd 2>/dev/null || echo inactive)"
  echo "adjtimex:   $(sudo /home/hpckey/bin/freqctl 2>/dev/null || echo 'freqctl 未编译')"
  echo "sshd:       $(sudo sshd -T 2>/dev/null | grep -i '^allowusers' || echo 无限制)"
  echo "orig 记录:  $([[ -f $REC_DIR/orig.env ]] && echo 有 || echo 无)"
}

restore() {
  if [[ ! -f "$REC_DIR/orig.env" ]]; then
    echo "没有原值记录 ($REC_DIR/orig.env)，拒绝恢复"; exit 1
  fi
  # shellcheck disable=SC1090
  source "$REC_DIR/orig.env"
  sudo cpupower -c all frequency-set -g "${governor:-performance}" || log "WARN: 恢复治理器失败"
  for f in /sys/devices/system/cpu/cpu[0-9]*/cpufreq/scaling_max_freq; do
    echo "${scaling_max_freq:-$(sysfs_get 0 scaling_max_freq)}" | sudo tee "$f" > /dev/null
  done
  for f in /sys/devices/system/cpu/cpu[0-9]*/cpufreq/scaling_min_freq; do
    echo "${scaling_min_freq:-$(sysfs_get 0 scaling_min_freq)}" | sudo tee "$f" > /dev/null
  done
  log "治理器与频率范围已恢复（${governor:-performance}）"
  if [[ "${chronyd_active:-inactive}" == "active" ]]; then
    sudo systemctl start chronyd && log "chronyd 已恢复"
  fi
  if [[ "${timesyncd_active:-inactive}" == "active" ]]; then
    sudo systemctl start systemd-timesyncd && log "systemd-timesyncd 已恢复"
  fi
  if [[ "${sshd_dropin_present:-no}" == "yes" ]]; then
    : # 原值本身就有 drop-in 才保留，否则删除
  fi
  sudo rm -f /etc/ssh/sshd_config.d/99-allow-hpckey.conf
  sudo systemctl reload ssh 2>/dev/null || sudo systemctl reload sshd 2>/dev/null || true
  log "sshd 限制已回滚（drop-in 删除）"
  log "restore 完成"
}

case "$ACTION" in
  setup)   setup ;;
  status)  status ;;
  restore) restore ;;
  *) echo "未知动作: $ACTION（setup|status|restore）"; exit 2 ;;
esac
exit 0