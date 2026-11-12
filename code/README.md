# Code guide

Run scripts from the package root. run_all.R executes steps 1-6 and the check in this order (step 3
only with RUN_BOOTSTRAP_WEIGHTS=TRUE); step 0 is run separately.

| Step | Script | Output |
|---|---|---|
| 0 | 00_prepare_data.R (optional; not run by run_all.R) | rebuilt data files in data/_prepared/, compared with data/ |
| 1 | 01_estimate_weights.R | results/point_estimates/SPEC.rds |
| 2 | 02_irf_bands.R | point IRFs and optional fresh bands |
| 3 | 03_bootstrap_weights.R (optional) | fresh bootstrap and expanding-window results |
| 4 | 04a_figures_irf.R, 04b_figures_series_weights.R | main figures |
| 5 | 05a_figures_irfboot.R, 05b_figures_bootstrap.R | appendix figures |
| 6 | 06_summary_numbers.R | results/summary_numbers.txt (numbers quoted in the text) |
| Check | funcs/validate_package.R | package checks (data, settings, archive dimensions, presence of all outputs) |

On a laptop (AMD Ryzen 7 PRO 5850U, 32 GB RAM) the default run_all.R takes 4-9 minutes and step 0
takes 5-10 minutes.

Reusable functions are in funcs/. Cluster launchers for step 3 are in slurm_helper/.
funcs/test_fast_var.R is a stand-alone test (not run by run_all.R) that checks that the fast
least-squares VAR refit used in the residual bootstrap matches vars::VAR to within 1e-9:
`Rscript code/funcs/test_fast_var.R`.

## Estimation settings

- VAR(12) with a constant and climate ordered first.
- Cholesky order: climate, INDPRO, PCEPI, UNRATE, SHADOWRATE, EXUSUKx, HOUST.
- Sample: January 1975-April 2025, with January 2020-June 2021 removed.
- Weights are non-negative and sum to one. NlcOptim::solnl starts from equal weights.
- Stored penalties come from five-fold contiguous block CV. Settings are in config.R.
- CV uses the calendar through May 2025. The final estimate ends in April 2025.
- A lambda is retained only when all five folds finish.

## Bootstrap and figures

- Residual bootstrap: block length 8, 1,000 outer draws, 100 inner draws, 200 standard-error
  refits, 68% percentile-t intervals and seed 123.
- Pairs bootstrap: block length 13, 1,000 draws, 68% basic intervals and L'Ecuyer-CMRG streams
  beginning at seed 123.
- Figure sizes, colours and fonts are set in funcs/theme.R and the figure functions. STIX Two Text is
  used when installed. Otherwise, R uses its default serif font.
