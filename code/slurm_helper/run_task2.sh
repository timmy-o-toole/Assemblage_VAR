#!/bin/bash
#SBATCH --job-name=avar_boot_t2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --output=slurm_task2_%j.log
#
# Task 2 -- expanding-window re-estimation of the point weights at all 216 origins (no bootstrap).
# Generic SLURM launcher. Site-specific environment setup can be supplied through
# R_ENV_SETUP; otherwise Rscript must already be available on PATH.
#
# Submit from the package root (the directory containing config.R):
#   sbatch code/slurm_helper/run_task2.sh                      # all six specifications, one after the other
#   SPECS="comp_noaa rank_noaa" sbatch code/slurm_helper/run_task2.sh
#   sbatch --array=0-5 code/slurm_helper/run_task2.sh          # one job per specification (parallel)
# Output: results/bootstrap/_recomputed/<spec>_task2.rds, logs in results/bootstrap/_logs/.
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
N_ORIGINS="${N_ORIGINS:-216}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  SPECS="${ALL_SPECS[$SLURM_ARRAY_TASK_ID]}"
else
  SPECS="${SPECS:-${ALL_SPECS[*]}}"
fi

mkdir -p results/bootstrap/_logs
for SPEC in $SPECS; do
  echo "============================================================"
  echo "TASK 2 | spec=${SPEC} | last ${N_ORIGINS} of 216 origins | $(date)"
  echo "Rscript code/03_bootstrap_weights.R ${SPEC} 2 0 ${N_ORIGINS}"
  echo "============================================================"
  Rscript code/03_bootstrap_weights.R "$SPEC" 2 0 "$N_ORIGINS" 2>&1 | tee "results/bootstrap/_logs/${SPEC}_task2.log"
done
echo "TASK 2 finished $(date)"
