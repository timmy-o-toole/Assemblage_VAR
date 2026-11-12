## =============================================================================
## 04a_figures_irf.R -- main-text impulse-response panels (Section 4).
##
## For every specification, draws the responses to a climate shock of the
## assembled index (red, fixed-weights residual-bootstrap 68% band) against the
## equal-weight headline ACI benchmark (blue, its own band).
##
## Reads  results/irf_bands/<spec>.rds and results/irf_bands/benchmark.rds
## Writes figures/irf_<spec>.png   (3 x 9.5 inch, 600 dpi, as the published files)
##
## Usage:  Rscript code/04a_figures_irf.R [spec_id ...]
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
  band_file <- require_result(irf_result_file(id), paste0("IRF band result for ", id))
  bands <- readRDS(band_file)
  if (is.null(bands$bands)) stop("IRF band result contains no bands: ", band_file, call. = FALSE)

  df <- irf_plot_data(bands, benchmark, spec_name = spec$label, bench_name = "Headline ACI")
  p  <- draw_irf_panel(df, y_scales = irf_y_scales(spec$rank),
                       labels = irf_panel_labels(spec$level, spec$rank))
  out <- file.path(PATHS$figures, paste0("irf_", id, ".png"))
  save_figure(p, out, width = 3, height = 9.5, units = "in", dpi = 600)
  cat("== ", id, ": saved ", file.path("figures", basename(out)), "\n", sep = "")
}
