## =============================================================================
## 05a_figures_irfboot.R -- appendix impulse-response panels with bands
## from the weight-reestimating pairs bootstrap.
##
## The band of the assembled index is recomputed from the VAR coefficient
## matrices and residuals saved in results/bootstrap/<spec>_task1.rds: the point
## response from var_B/var_resid, one response per converged draw from
## boot$var_B/boot$var_resid (Cholesky, climate index first, horizon 36,
## cumulated except SHADOWRATE), then 68% basic (reverse-percentile) intervals.
## The plotted point response is the fixed-weights
## one of results/irf_bands/<spec>.rds; only its lower/upper bands are replaced.
## The benchmark (blue) keeps its residual-bootstrap band.
##
## Reads  the archived results/bootstrap/<spec>_task1.rds (or, with USE_RECOMPUTED=TRUE,
##        results/bootstrap/_recomputed/<spec>_task1.rds), results/irf_bands/<spec>.rds,
##        results/irf_bands/benchmark.rds
## Writes figures/irfboot_<spec>.png   (3 x 9.5 inch, 600 dpi)
##
## Usage:  Rscript code/05a_figures_irfboot.R [spec_id ...]
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "theme.R"))
source(file.path(ROOT, "code", "funcs", "fig_irf.R"))

args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args)) args else names(SPECS)
unknown <- setdiff(ids, names(SPECS))
if (length(unknown)) stop("unknown specification id(s): ", paste(unknown, collapse = ", "))

benchmark_file <- require_result(irf_result_file("benchmark"), "benchmark IRF band result")
benchmark <- readRDS(benchmark_file)

for (id in ids) {
  spec <- SPECS[[id]]
  task_file <- boot_result_file(id, 1)
  band_file <- irf_result_file(id)
  require_result(task_file, paste0("task-1 bootstrap result for ", id))
  require_result(band_file, paste0("IRF band result for ", id))
  bands <- readRDS(band_file)
  if (is.null(bands$bands)) stop("IRF band result contains no bands: ", band_file, call. = FALSE)
  task1 <- readRDS(task_file)
  cat("== ", id, ": appendix inputs ", boot_result_source(task_file), "\n", sep = "")

  ## responses from the saved VARs: point estimate and converged bootstrap draws
  cumulate <- vapply(bands$bands, function(b) isTRUE(b$cumulated), logical(1))
  pb <- pairs_bootstrap_irf(task1, var_order = bands$var_order, p = CFG$lags, horizon = CFG$irf_horizon,
                            cumulate = cumulate, macro_order = CFG$macro_order)
  cat(sprintf("== %s: %d of %d usable bootstrap draws\n", id, pb$n_draws, pb$n_total))

  ## 68% reverse-percentile bands around the recomputed point response, plotted with the
  ## fixed-weights point response used in the appendix panel
  ci <- pairs_bootstrap_bands(pb$point, pb$draws, conf = CFG$ci_level, reverse_percentile = TRUE)
  bands_boot <- replace_bands(bands, ci$lower, ci$upper, source = "weight-reestimating pairs bootstrap")

  df <- irf_plot_data(bands_boot, benchmark, spec_name = spec$label, bench_name = "Headline ACI")
  p  <- draw_irf_panel(df, y_scales = irf_y_scales(spec$rank, appendix = TRUE),
                       labels = irf_panel_labels(spec$level, spec$rank, appendix = TRUE))
  out <- file.path(PATHS$figures, paste0("irfboot_", id, ".png"))
  save_figure(p, out, width = 3, height = 9.5, units = "in", dpi = 600)
  cat("   saved ", file.path("figures", basename(out)), "\n", sep = "")
}
