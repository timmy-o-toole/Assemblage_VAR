## =============================================================================
## 04b_figures_series_weights.R -- assembled-index series and weight
## figures of the main text (Figures 3-4).
##
## Reads  results/point_estimates/<spec>.rds (weights), the seasonally adjusted level
##        panels data/aci_national_level.csv, aci_regional_level.csv, noaa_level.csv and
##        the transformed panels (rank order) through code/funcs/load_data.R and
##        code/funcs/appendix_summaries.R, data/labels_*.csv, data/nber_recessions.csv
## Writes figures/fitted_aci.png        headline ACI vs national / regional
##                                      component index and regional rank index
##        figures/fitted_noaa.png       headline ACI vs national / NOAA component
##                                      index and NOAA rank index
##        figures/weights_comp_aci.png  national and regional component weights
##        figures/weights_rank_aci.png  regional rank profile (8-rank trapezoid)
##        figures/weights_comp_noaa.png NOAA weights by region x variable cell
##        figures/weights_rank_noaa.png NOAA rank profile (2-rank average)
## Sizes and dpi are those of the published files (46 x 18 cm at 150 dpi for the
## series; 300 x 205 mm at 150 dpi, scale 1.1, for the weights).
##
## Usage:  Rscript code/04b_figures_series_weights.R
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "theme.R"))
source(file.path(ROOT, "code", "funcs", "load_data.R"))
source(file.path(ROOT, "code", "funcs", "appendix_summaries.R"))   # smoothed level panels (as Appendix A.2)
source(file.path(ROOT, "code", "funcs", "fig_series_weights.R"))

read_weights <- function(id) readRDS(file.path(PATHS$estimates, paste0(id, ".rds")))$weights

weights <- lapply(stats::setNames(nm = names(SPECS)), read_weights)
weights <- lapply(weights, function(w) {
  if (any(w < -1e-10)) stop("materially negative estimated weight found", call. = FALSE)
  w[w < 0] <- 0
  w / sum(w)
})

labels_regional <- utils::read.csv(file.path(PATHS$data, "labels_aci_regional.csv"), stringsAsFactors = FALSE, check.names = FALSE)
labels_noaa     <- utils::read.csv(file.path(PATHS$data, "labels_noaa.csv"),         stringsAsFactors = FALSE, check.names = FALSE)

## ---- assembled indices --------------------------------------------------------------------------
## Five-year (60-month) moving averages of the standardised level panels, labelled by the last
## month of the window (first value 1965-12); see smoothed_level_panel(). Rank panels are the
## smoothed levels sorted by the rank order of the transformed panel at the same date.
level_panel <- c("national", "regional", "noaa")
idx <- function(id) {
  panel <- level_panel[SPECS[[id]]$level]
  X <- if (SPECS[[id]]$rank) smoothed_rank_panel(panel, PATHS$data) else smoothed_component_panel(panel, PATHS$data)
  assembled_index(X, weights[[id]])
}
smooth <- function(x, k) { out <- trailing_mean(x, k); names(out) <- names(x); out }

headline <- smoothed_aci_headline(PATHS$data)       # equal-weight ACI, all six national components
comp_nat <- idx("comp_national")

## ACI figure: regional component index, regional rank index (extra 24-month average), national,
## headline. Window: the full ACI sample (1961-01 to 2025-05) as in the paper, no recession
## shading; the smoothed series start in 1965-12 because the level panels start in 1961-01.
series_aci <- list(
  "Assemblage (comps 2)" = idx("comp_regional"),
  "Assemblage (ranks 2)" = smooth(idx("rank_regional"), 24L),
  "Assemblage (comps 1)" = comp_nat,
  "ACI"                  = headline
)
colours_aci <- c("Assemblage (comps 2)" = "#ef3b2c", "Assemblage (ranks 2)" = "midnightblue",
                 "Assemblage (comps 1)" = "firebrick", "ACI" = "#1B9E77")
labels_aci  <- expression(MRCI[scriptstyle(comps)]^{scriptstyle(Region)}, MRCI[scriptstyle(ranks)]^{scriptstyle(Region)},
                          MRCI[scriptstyle(comps)]^{scriptstyle(Nation)}, ACI)
window_aci <- rownames(read_panel(file.path(PATHS$data, "aci_national_level.csv")))
p <- draw_series_plot(series_plot_frame(series_aci, window_aci), colours_aci, labels_aci)
save_figure(p, file.path(PATHS$figures, "fitted_aci.png"), width = 46, height = 18, units = "cm", dpi = 150)
cat("== saved figures/fitted_aci.png\n")

## NOAA figure: NOAA component index, NOAA rank index (extra 12-month average), national,
## headline. Window: the months of the smoothed NOAA panel within the ACI sample (1965-12 to
## 2025-05), with NBER recession shading.
series_noaa <- list(
  "Assemblage (comps 3)" = idx("comp_noaa"),
  "Assemblage (ranks 3)" = smooth(idx("rank_noaa"), 12L),
  "Assemblage (comps 1)" = comp_nat,
  "ACI"                  = headline
)
colours_noaa <- c("Assemblage (comps 3)" = "#fc9272", "Assemblage (ranks 3)" = "#2171b5",
                  "Assemblage (comps 1)" = "firebrick", "ACI" = "#1B9E77")
labels_noaa_expr <- expression(MRCI[scriptstyle(comps)]^{scriptstyle(NOAA)}, MRCI[scriptstyle(ranks)]^{scriptstyle(NOAA)},
                               MRCI[scriptstyle(comps)]^{scriptstyle(Nation)}, ACI)
window_noaa <- intersect(names(series_noaa[[1]]), names(headline))
p <- draw_series_plot(series_plot_frame(series_noaa, window_noaa), colours_noaa, labels_noaa_expr,
                      recession_frame(PATHS$data, window_noaa[1]))
save_figure(p, file.path(PATHS$figures, "fitted_noaa.png"), width = 46, height = 18, units = "cm", dpi = 150)
cat("== saved figures/fitted_noaa.png\n")

## ---- weights ------------------------------------------------------------------------------------
p <- draw_weights_comp_aci(weights$comp_regional, weights$comp_national, labels_regional)
save_figure(p, file.path(PATHS$figures, "weights_comp_aci.png"), width = 300, height = 205, units = "mm", dpi = 150, scale = 1.1)
cat("== saved figures/weights_comp_aci.png\n")

p <- draw_weights_rank(smooth_rank_weights_trapezoid(weights$rank_regional),
                       legend_label = expression(MRCI[scriptstyle(ranks)]^{scriptstyle(Region)}),
                       text_size = 36.4, aspect_ratio = 0.572)
save_figure(p, file.path(PATHS$figures, "weights_rank_aci.png"), width = 300, height = 205, units = "mm", dpi = 150, scale = 1.1)
cat("== saved figures/weights_rank_aci.png\n")

p <- draw_weights_comp_noaa(weights$comp_noaa, labels_noaa)
save_figure(p, file.path(PATHS$figures, "weights_comp_noaa.png"), width = 300, height = 205, units = "mm", dpi = 150)
cat("== saved figures/weights_comp_noaa.png\n")

p <- draw_weights_rank(smooth_rank_weights_pair(weights$rank_noaa),
                       legend_label = expression(MRCI[scriptstyle(ranks)]^{scriptstyle(NOAA)}),
                       text_size = 35.6, aspect_ratio = 0.575)
save_figure(p, file.path(PATHS$figures, "weights_rank_noaa.png"), width = 300, height = 205, units = "mm", dpi = 150, scale = 1.1)
cat("== saved figures/weights_rank_noaa.png\n")
