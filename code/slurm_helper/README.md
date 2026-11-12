# Step 3 on SLURM

03_bootstrap_weights.R runs two optional jobs for each specification:

| Task | Work |
|---|---|
| 1 | Final weights and 1,000 pairs-bootstrap draws |
| 2 | Weights at 216 expanding-window origins |

The launchers request 2 CPUs, 4 GB per CPU and 120 hours. Requirements: the complete package,
Rscript on PATH, and the NlcOptim and zoo packages. Use R_ENV_SETUP for site-specific modules.

## Submit

    cd <package root>
    sbatch code/slurm_helper/run_task1.sh
    sbatch code/slurm_helper/run_task2.sh

Optional subsets and arrays:

    SPECS="comp_noaa rank_noaa" sbatch code/slurm_helper/run_task1.sh
    R_ENV_SETUP="/path/to/r-setup.sh" sbatch code/slurm_helper/run_task1.sh
    sbatch --array=0-5 code/slurm_helper/run_task1.sh

Outputs are written to results/bootstrap/_recomputed/. Package archives in results/bootstrap/ are
never overwritten. Set USE_RECOMPUTED=TRUE when drawing figures from fresh results.

## Smoke runs

    Rscript code/03_bootstrap_weights.R comp_national 1 5
    Rscript code/03_bootstrap_weights.R comp_national 2 0 3

Reduced runs go to results/bootstrap/_smoke/. Random streams use L'Ecuyer-CMRG and seed 123. Origin
o uses seed 123 + o. Optimizer output may vary with the numerical platform.
