#!/bin/bash
#SBATCH --job-name=avar_boot_t1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --output=slurm_task1_%j.log
#
# Task 1 -- weight-reestimating pairs bootstrap at the final origin, B = 1000 draws.
# Generic SLURM launcher. Site-specific environment setup can be supplied through
# R_ENV_SETUP; otherwise Rscript must already be available on PATH.
#
# Submit from the package root (the directory containing config.R):
#   sbatch code/slurm_helper/run_task1.sh                      # all six specifications, one after the other
#   SPECS="comp_national rank_national" sbatch code/slurm_helper/run_task1.sh
#   sbatch --array=0-5 code/slurm_helper/run_task1.sh          # one job per specification (parallel)
# Output: results/bootstrap/_recomputed/<spec>_task1.rds, logs in results/bootstrap/_logs/.
# The archived outputs in results/bootstrap/ are never overwritten.

set -euo pipefail

if [[ -n "${R_ENV_SETUP:-}" ]]; then
  test -f "$R_ENV_SETUP" || { echo "R_ENV_SETUP not found: $R_ENV_SETUP"; exit 1; }
  # shellcheck source=/dev/null
  source "$R_ENV_SETUP"
fi
command -v Rscript >/dev/null 2>&1 || { echo "Rscript is not available on PATH"; exit 1; }

ROOT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$ROOT_DIR"
test -f config.R || { echo "config.R not found in $ROOT_DIR -- submit from the package root"; exit 1; }

ALL_SPECS=(comp_national comp_regional comp_noaa rank_national rank_regional rank_noaa)
B="${B:-1000}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  SPECS="${ALL_SPECS[$SLURM_ARRAY_TASK_ID]}"
else
  SPECS="${SPECS:-${ALL_SPECS[*]}}"
fi

mkdir -p results/bootstrap/_logs
for SPEC in $SPECS; do
  echo "============================================================"
  echo "TASK 1 | spec=${SPEC} | B=${B} | origin 216 (full sample) | $(date)"
  echo "Rscript code/03_bootstrap_weights.R ${SPEC} 1 ${B}"
  echo "============================================================"
  Rscript code/03_bootstrap_weights.R "$SPEC" 1 "$B" 2>&1 | tee "results/bootstrap/_logs/${SPEC}_task1.log"
done
echo "TASK 1 finished $(date)"
