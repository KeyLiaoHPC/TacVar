#!/usr/bin/env bash
# Run TacVar tests on ARM host c920bn1 via ssh (shared NFS workspace).
#   bash src/measure/tests/run_arm_tests.sh           # unit + suite smoke (<=4 cores)
#   bash src/measure/tests/run_arm_tests.sh --tf-c920 # CG.C in-situ TF, np=128
set -euo pipefail
HOST="${TACVAR_ARM_HOST:-c920bn1}"
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
REMOTE_ROOT="$ROOT"
MODE="${1:-}"

ssh -o BatchMode=yes -o ConnectTimeout=10 "$HOST" "test -d '$REMOTE_ROOT'" \
  || { echo "ERROR: $HOST cannot see $REMOTE_ROOT"; exit 1; }

if [[ "$MODE" == "--tf-c920" ]]; then
  ssh "$HOST" bash -s <<EOF
set -euo pipefail
export PATH=/home/hpckey/01-App/openmpi-5.0.8/bin:\$PATH
export LD_LIBRARY_PATH=/home/hpckey/01-App/openmpi-5.0.8/lib:/home/hpckey/01-App/papi-7.2.0b2/lib:\${LD_LIBRARY_PATH:-}
export PAPI_HOME=/home/hpckey/01-App/papi-7.2.0b2
export MPI_HOME=/home/hpckey/01-App/openmpi-5.0.8
export TACVAR_NSTP_ARM=10
cd '$REMOTE_ROOT'
uname -m
bash suites/NPB3.4.4/NPB3.4-MPI/scripts/test_tf_c920bn1.sh
EOF
  exit 0
fi

ssh "$HOST" bash -s <<EOF
set -euo pipefail
export PATH=/home/hpckey/01-App/openmpi-5.0.8/bin:\$PATH
export LD_LIBRARY_PATH=/home/hpckey/01-App/openmpi-5.0.8/lib:/home/hpckey/01-App/papi-7.2.0b2/lib:\${LD_LIBRARY_PATH:-}
export PAPI_HOME=/home/hpckey/01-App/papi-7.2.0b2
export MPI_HOME=/home/hpckey/01-App/openmpi-5.0.8
export TACVAR_NSTP_ARM=10
cd '$REMOTE_ROOT'
uname -m
ls /sys/module/ph_enable_pmu >/dev/null 2>&1 && echo kmod=yes || echo kmod=no
bash src/measure/tests/run_unit_tests.sh
bash src/measure/tests/run_backend_smoke.sh
# Small suite smokes (<=4 cores); use arch-local build dirs via force rebuild in scripts
unset TACVAR_DATA_DIR || true
bash suites/NPB3.4.4/test_tacvar.sh --run-smoke
echo "arm: NPB smoke finished"
unset TACVAR_DATA_DIR || true
bash suites/lmbench/scripts/test_tacvar.sh --run-smoke
echo "arm: lmbench smoke finished"
EOF
