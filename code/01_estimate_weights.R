## =============================================================================
## 01_estimate_weights.R -- Step 1: aggregation weights for the six
## specifications (Section 3 of the paper).
##
## For each specification: load the panel, (optionally) select the penalty by
## block cross-validation, solve the constrained problem, and save the weights,
## the assembled climate index and the VAR data to results/point_estimates/.
##
## Usage:  Rscript code/01_estimate_weights.R [spec_id ...]
##         RUN_CV=TRUE Rscript code/01_estimate_weights.R   (re-select lambda; slow)
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "load_data.R"))
source(file.path(ROOT, "code", "funcs", "estimation.R"))
source(file.path(ROOT, "code", "funcs", "cv.R"))

args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args)) args else names(SPECS)
set.seed(SEED)

for (id in ids) {
  spec <- SPECS[[id]]
  cat(sprintf("\n== %s (%s) ==\n", id, spec$label))
  t0 <- Sys.time()
  cfg <- CFG; cfg$climate_name <- "climate"
  d   <- load_specification(spec, PATHS$data, cfg)
  est <- prepare_estimation(d, CFG$lags)
  cat(sprintf("  sample %s..%s, %d months after removing %d covid months; %d lag-stack rows, %d gap-valid pairs; K = %d\n",
              min(d$dates), max(d$dates), length(d$dates), length(cfg$covid_dates),
              est$n_use, length(est$valid_pairs_upto(est$n_use)), est$K))

  cv <- NULL
  lambda <- spec$lambda
  if (RUN_CV) {
    cv <- cross_validate(est, d, spec, EST$lambda_grid(spec$rank), EST, n_cores = N_CORES)
    lambda <- cv$lambda
    cat(sprintf("  CV selected lambda = %.6g (grid index %d; published %.6g, index %d)\n",
                lambda, cv$grid_index, spec$lambda, spec$grid_index))
  } else {
    cat(sprintf("  using the published cross-validated lambda = %.6g\n", lambda))
  }

  fit <- solve_weights(est, d, spec, lambda, EST, mode = "final")
  w <- fit$weights
  cat(sprintf("  solved in %d criterion evaluations, %.1f s; sum(w) = %.6f, %d of %d weights > 1e-6\n",
              fit$n_eval, as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(w), sum(w > 1e-6), length(w)))
  top <- sort(w, decreasing = TRUE)[1:min(5, length(w))]
  cat("  largest weights:", paste(sprintf("%s=%.3f", names(top), top), collapse = ", "), "\n")

  out <- list(
    id = id, spec = spec, lambda = lambda, lambda_source = if (RUN_CV) "cross-validation (this run)" else "published cross-validated value",
    weights = w, objective = fit$objective,
    index = assembled_index(d, w),              # w'X_t over the estimation sample (covid months removed)
    benchmark = d$bench[[d$climate_name]],       # equal-weight (headline) ACI on the same dates
    dates = d$dates, dates_full = d$dates_full, gap_dummy = d$gap_dummy,
    macro = d$macro, bench = d$bench, comps = d$comps, comps_perm = d$comps_perm,
    var_order = d$var_order, lags = CFG$lags, end_test_date = d$end_test_date,
    cv = if (!is.null(cv)) cv[c("lambda", "grid_index", "fold_mean", "grid", "fold_sizes", "n_pool")] else NULL,
    session = list(time = Sys.time(), R = R.version.string, NlcOptim = as.character(packageVersion("NlcOptim")))
  )
  saveRDS(out, file.path(PATHS$estimates, paste0(id, ".rds")))
  cat("  saved", file.path("results/point_estimates", paste0(id, ".rds")), "\n")
}
