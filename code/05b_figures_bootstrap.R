## Appendix A.1/A.2/B figures from supplied or freshly recomputed stability results.
## This script performs no estimation.  Numerical preparation lives in
## code/funcs/appendix_summaries.R so it can be tested independently of rendering.
source("config.R")
source(file.path(ROOT, "code", "funcs", "load_data.R"))
source(file.path(ROOT, "code", "funcs", "theme.R"))
source(file.path(ROOT, "code", "funcs", "appendix_summaries.R"))
source(file.path(ROOT, "code", "funcs", "appendix_style.R"))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

lab_reg <- read.csv(file.path(PATHS$data, "labels_aci_regional.csv"), stringsAsFactors = FALSE)
lab_noaa <- read.csv(file.path(PATHS$data, "labels_noaa.csv"), stringsAsFactors = FALSE)
## --- assembled-index panels (Appendix A.2): colours, labels and shading of the published figure
index_colours <- c(comp_national = "firebrick", comp_regional = "#ef3b2c", comp_noaa = "#fc4e2a",
                   rank_national = "#08519c", rank_regional = "midnightblue", rank_noaa = "#2171b5")
index_labels <- list(
  comp_national = expression(MRCI[scriptstyle(comps)]^{scriptstyle(Nation)}),
  comp_regional = expression(MRCI[scriptstyle(comps)]^{scriptstyle(Region)}),
  comp_noaa     = expression(MRCI[scriptstyle(comps)]^{scriptstyle(NOAA)}),
  rank_national = expression(MRCI[scriptstyle(ranks)]^{scriptstyle(Nation)}),
  rank_regional = expression(MRCI[scriptstyle(ranks)]^{scriptstyle(Region)}),
  rank_noaa     = expression(MRCI[scriptstyle(ranks)]^{scriptstyle(NOAA)}))

lighten <- function(col, amount = 0.55) {
  v <- grDevices::col2rgb(col) / 255
  v <- v + (1 - v) * amount
  grDevices::rgb(v[1], v[2], v[3])
}

## NBER recessions overlapping a plotting window, with the first one clipped to the window start.
nber_bands <- function(window) {
  rec <- read.csv(file.path(PATHS$data, "nber_recessions.csv"), stringsAsFactors = FALSE)
  rec$start <- as.Date(rec$start); rec$end <- as.Date(rec$end)
  open <- rec$start <= window[1] & rec$end >= window[1]
  rec$start[open] <- window[1]
  rec[rec$end >= window[1] & rec$start <= window[2], , drop = FALSE]
}

rank_groups <- function(K) {
  paste0("Q", pmin(5L, ceiling(seq_len(K) * 5 / K)))
}

for (id in names(SPECS)) {
  t1f <- boot_result_file(id, 1)
  t2f <- boot_result_file(id, 2)
  pef <- file.path(PATHS$estimates, paste0(id, ".rds"))
  require_result(t1f, paste0("task-1 bootstrap result for ", id))
  require_result(t2f, paste0("task-2 expanding-window result for ", id))
  require_result(pef, paste0("point estimate for ", id))
  point <- readRDS(pef); task1 <- readRDS(t1f); task2 <- readRDS(t2f)
  cat("== ", id, ": appendix inputs ", boot_result_source(t1f), "\n", sep = "")
  ## The plotted point is the main-text estimate; the intervals come from the bootstrap draws.
  names(task1$weights) <- names(point$weights)
  cfg <- CFG; cfg$climate_name <- "climate"
  d <- load_specification(SPECS[[id]], PATHS$data, cfg)

  ## Bootstrap weights: family intervals in component space, fan chart in rank space.
  ws <- bootstrap_weight_statistics(id, point$weights, task1, lab_reg, lab_noaa)
  if (startsWith(id, "comp_")) {
    boot_levels <- appendix_display_levels(ws$group, appendix_component_levels(id, "bootstrap"))
    ws$group <- factor(ws$group, levels = boot_levels)
    p <- ggplot(ws, aes(group, 100 * point, fill = group)) +
      geom_col(width = .72, alpha = .85) +
      geom_linerange(aes(ymin = 100 * q05, ymax = 100 * q95), linewidth = .7) +
      geom_linerange(aes(ymin = 100 * q16, ymax = 100 * q84), linewidth = 2.5) +
      geom_point(aes(y = 100 * q50), shape = 21, fill = "white", size = 2.2) +
      scale_fill_manual(values = APPENDIX_FAMILY_COLOURS[boot_levels], guide = "none") +
      labs(x = NULL, y = "Weights (in %)") + theme_publication(18)
  } else {
    ws$x <- seq_len(nrow(ws))
    p <- ggplot(ws, aes(x)) +
      geom_ribbon(aes(ymin = 100 * q05, ymax = 100 * q95), fill = "#C2D8F2") +
      geom_ribbon(aes(ymin = 100 * q16, ymax = 100 * q84), fill = "#3C68BB") +
      geom_line(aes(y = 100 * point), colour = "#001330", linewidth = 1) +
      labs(x = "Rank", y = "Weights (in %)") + theme_publication(18)
  }
  save_figure(p, file.path(PATHS$figures, paste0("boot_weights_", id, ".png")),
              width = 300, height = 205, units = "mm", dpi = 150, scale = 1.1)

  ## Bootstrap assembled-index envelope (Appendix A.2).  Drawn on five-year moving averages of the
  ## seasonally adjusted level panels, as in the published figure; see bootstrap_index_paths().
  bi <- bootstrap_index_paths(id, SPECS[[id]], PATHS$data, point$weights, task1)
  main <- unname(index_colours[id]); light <- lighten(main, 0.55)
  rec <- nber_bands(range(bi$date))
  line_labels <- list(index_labels[[id]], "ACI"); names(line_labels) <- c("point", "aci")
  p <- ggplot(bi, aes(date)) +
    geom_rect(data = rec, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "blue2", alpha = 0.05) +
    geom_ribbon(aes(ymin = q05, ymax = q95, fill = "90% CI"), alpha = 0.55) +
    geom_ribbon(aes(ymin = q16, ymax = q84, fill = "68% CI"), alpha = 0.55) +
    geom_hline(yintercept = 0, linewidth = 0.6, colour = "grey29") +
    geom_line(aes(y = benchmark, colour = "aci"), linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = point, colour = "point"), linewidth = 1.15, na.rm = TRUE) +
    scale_colour_manual(name = NULL, breaks = c("point", "aci"),
                        values = c(point = main, aci = "#1B9E77"), labels = line_labels) +
    scale_fill_manual(name = NULL, breaks = c("68% CI", "90% CI"),
                      values = c(`68% CI` = main, `90% CI` = light)) +
    scale_x_date(breaks = scales::pretty_breaks(n = 8), date_labels = "%Y",
                 expand = c(0.01, 0.01)) +
    scale_y_continuous(labels = function(x) sprintf("% .2f", x),
                       expand = expansion(mult = c(0.04, 0.09))) +
    labs(x = NULL, y = NULL) + theme_publication(22) +
    theme(legend.position = "inside", legend.position.inside = c(0.012, 0.95),
          legend.justification = c(0, 1),
          legend.direction = "horizontal", legend.box = "vertical",
          legend.spacing.y = unit(0.02, "cm"), legend.margin = margin(2, 6, 2, 4),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.key = element_rect(fill = "transparent", colour = NA))
  save_figure(p, file.path(PATHS$figures, paste0("boot_index_", id, ".png")),
              width = 46, height = 18.5, units = "cm", dpi = 150)

  ## Expanding-window profiles.  COVID months are inserted at December-2019
  ## values before applying the paper's centred partial 12-month mean.
  ew <- expanding_weight_calendar(task2$W_mat, smooth = TRUE)
  if (startsWith(id, "comp_")) {
    groups <- appendix_groups(id, point$weights, lab_reg, lab_noaa)
    lev <- appendix_display_levels(groups, appendix_component_levels(id, "expanding"))
    sh <- aggregate_weight_rows(ew, groups, lev)
  } else {
    ## Quintiles of the displayed (rank-smoothed) profile, as in the paper.
    ew <- expanding_weight_calendar(t(apply(task2$W_mat, 1L, smooth_rank_profile, id = id)), smooth = TRUE)
    groups <- rank_groups(ncol(ew)); lev <- APPENDIX_RANK_ORDER
    sh <- aggregate_weight_rows(ew, groups, lev)
  }
  sh <- sh / rowSums(sh)
  long <- data.frame(date = rep(as.Date(rownames(sh)), times = ncol(sh)),
                     group = factor(rep(colnames(sh), each = nrow(sh)), levels = colnames(sh)),
                     weight = 100 * as.numeric(sh))
  cols <- if (startsWith(id, "comp_")) APPENDIX_FAMILY_COLOURS[lev] else APPENDIX_RANK_COLOURS[lev]
  p <- ggplot(long, aes(date, weight, fill = group)) +
    geom_area(position = position_stack(reverse = TRUE), alpha = APPENDIX_EXPANDING_ALPHA) +
    scale_fill_manual(values = cols, breaks = lev) +
    labs(x = NULL, y = "Weights (in %)", fill = NULL) + theme_publication(18)
  save_figure(p, file.path(PATHS$figures, paste0("expanding_weights_", id, ".png")),
              width = 300, height = 205, units = "mm", dpi = 150, scale = 1.1)
  cat("== ", id, ": bootstrap-weight, index, and expanding-window figures saved\n", sep = "")
}

## Appendix B: rank-implied NOAA usage by region-family cell.
rank_noaa_pe <- file.path(PATHS$estimates, "rank_noaa.rds")
comp_noaa_pe <- file.path(PATHS$estimates, "comp_noaa.rds")
rank_noaa_t1 <- boot_result_file("rank_noaa", 1)
comp_noaa_t1 <- boot_result_file("comp_noaa", 1)
required_usage <- c(rank_noaa_pe, comp_noaa_pe, rank_noaa_t1, comp_noaa_t1)
if (all(file.exists(required_usage))) {
  cat("== rank_usage_noaa: appendix inputs ", boot_result_source(rank_noaa_t1), "\n", sep = "")
  rr <- readRDS(rank_noaa_pe)
  rc <- readRDS(comp_noaa_pe)
  u <- rank_usage_noaa(rr$weights, rc$weights, rr$comps_perm, lab_noaa)$cells
  u$family <- factor(u$family, levels = APPENDIX_NOAA_BOOT_ORDER)
  u$region <- factor(u$region, levels = APPENDIX_REGION_ORDER)
  usage_keys <- paste(as.character(u$family), as.character(u$region), sep = "::")
  u$key <- factor(usage_keys, levels = appendix_usage_levels(usage_keys))
  u$x <- as.integer(u$key)
  family_positions <- split(u$x, u$family, drop = TRUE)
  family_mid <- vapply(family_positions, function(x) mean(range(x)), numeric(1))
  family_end <- vapply(family_positions, max, integer(1))
  p <- ggplot(u, aes(x, 100 * usage, fill = region)) +
    geom_col(width = .8, alpha = APPENDIX_USAGE_ALPHA) + geom_hline(yintercept = 0) +
    scale_x_continuous(breaks = unname(family_mid), labels = names(family_mid),
                       minor_breaks = head(unname(family_end), -1) + .5,
                       expand = expansion(add = c(.5, .5))) +
    scale_fill_manual(values = APPENDIX_REGION_COLOURS, breaks = APPENDIX_REGION_ORDER) +
    labs(x = NULL, y = "Deviation from equal weight (in %)", fill = NULL) +
    theme_publication(16) +
    theme(axis.text.x = element_text(angle = 0, margin = margin(t = 6)),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(colour = "gray60", linetype = 2),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(12, 0, 0, 0))
  save_figure(p, file.path(PATHS$figures, "rank_usage_noaa.png"),
              width = 300, height = 185, units = "mm", dpi = 140)
  cat("== rank_usage_noaa: saved\n")
} else {
  stop("rank-usage input missing: ",
       paste(sub(ROOT, ".", required_usage[!file.exists(required_usage)], fixed = TRUE), collapse = ", "),
       call. = FALSE)
}
