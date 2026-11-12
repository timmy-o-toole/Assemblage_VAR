## =============================================================================
## run_all.R -- reproduce the package figures and numerical outputs.
##
##   Rscript run_all.R                       fast mode (default): uses the published cross-validated
##                                           penalties and the shipped impulse-response band objects;
##                                           re-estimates the weights, recomputes all point responses,
##                                           redraws every figure and writes the numbers quoted in
##                                           the paper to results/summary_numbers.txt. Allow about
##                                           10 minutes.
##   RUN_CV=TRUE Rscript run_all.R           also re-selects the penalties by block cross-validation
##                                           (30 x 5 constrained fits per specification; about 2 min
##                                           for the 5-series specifications up to 45-75 min for the
##                                           82-series NOAA specifications).
##   RUN_IRF_BOOTSTRAP=TRUE Rscript run_all.R   also recomputes the fixed-weights residual bootstrap
##                                           bands (1000 x 100 + 200 VAR refits per VAR; about 1 h
##                                           per VAR on one core, seven VARs).
##   RUN_BOOTSTRAP_WEIGHTS=TRUE Rscript run_all.R  also re-runs the appendix pairs bootstrap and the
##                                           expanding-window re-estimation for all six
##                                           specifications (cluster-scale; see the warning below).
##   USE_RECOMPUTED=TRUE ...                 draw all bootstrap-based figures and numbers only from
##                                           the corresponding _recomputed directories.
##   N_CORES=8 ...                           workers for the cross-validation.
##
## By default NOTHING expensive is recomputed: the penalties are the published cross-validated
## values, the residual-bootstrap bands are the shipped objects, and the appendix figures read the
## archived pairs-bootstrap and expanding-window estimation outputs in results/bootstrap/.
## Supplied results are never overwritten. Fresh IRF bands and appendix stability results are
## written to their respective _recomputed directories. Missing requested inputs are fatal.
## The package check (code/funcs/validate_package.R) runs last: it verifies the data, the settings,
## the supplied archives and the regenerated estimates and figures.
## =============================================================================
t_start <- Sys.time()
source("config.R")
run <- function(script, ...) {
  cat(sprintf("\n############ %s ############\n", script))
  if (!file.exists(file.path("code", script))) stop("required script not found: code/", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c(file.path("code", script), ...))
  if (status != 0) stop("step failed: ", script)
}
run("01_estimate_weights.R")
run("02_irf_bands.R")
if (RUN_BOOTSTRAP_WEIGHTS) {
  cat("\n!! RUN_BOOTSTRAP_WEIGHTS = TRUE: re-running 12 cluster-scale jobs (hours to days on one machine).\n",
      "!! Output goes to results/bootstrap/_recomputed/; the supplied outputs are not touched.\n",
      "!! On a cluster use code/slurm_helper/ instead. Add USE_RECOMPUTED=TRUE to plot the fresh run.\n", sep = "")
  for (id in names(SPECS)) {
    run("03_bootstrap_weights.R", id, "1", "1000")
    run("03_bootstrap_weights.R", id, "2", "0", "216")
  }
}
run("04a_figures_irf.R")
run("04b_figures_series_weights.R")
run("05a_figures_irfboot.R")
run("05b_figures_bootstrap.R")
run("06_summary_numbers.R")
run(file.path("funcs", "validate_package.R"))
cat(sprintf(paste0("\nrun_all.R finished in %.1f minutes (RUN_CV = %s, RUN_IRF_BOOTSTRAP = %s,\n",
                   "RUN_BOOTSTRAP_WEIGHTS = %s, USE_RECOMPUTED = %s; appendix inputs: %s)\n"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins")), RUN_CV, RUN_IRF_BOOTSTRAP,
            RUN_BOOTSTRAP_WEIGHTS, USE_RECOMPUTED,
            result_source(boot_result_file("comp_national", 1))))
