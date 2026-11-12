## =============================================================================
## code/funcs/fig_series_weights.R -- assembled-index time series and weight figures of
## the main text (Figures 3-4).
##
## Series figures and weight figures (component and rank branches), drawn from the
## package result objects with the layout of the published figures.
##
## Series: the assembled index is the point-estimate weight vector applied to
## the five-year (60-month) moving average of the standardised, seasonally
## adjusted LEVEL panel (smoothed_component_panel / smoothed_rank_panel in
## code/funcs/appendix_summaries.R, the construction of the published ACI and of
## the Appendix A.2 index figures), not to the monthly transformed panel the
## estimator uses. The rank index additionally carries the sliding average used
## in the paper (24 months for the ACI regional ranks, 12 months for the NOAA
## ranks).
##
## Requires: ggplot2, zoo, scales, ggthemes (through code/funcs/theme.R).
## =============================================================================
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

## ---- helpers ------------------------------------------------------------------------------

## Trailing moving average of length k; the first k-1 values are NA (as c(rep(NA, k-1),
## zoo::rollmean(x, k))).
trailing_mean <- function(x, k) {
  x <- as.numeric(x)
  if (k <= 1L) return(x)
  c(rep(NA_real_, k - 1L), zoo::rollmean(x, k))
}

## Assembled index w'X_t over a (untrimmed) component panel; named by date.
assembled_index <- function(comps, weights) {
  m <- as.matrix(comps)
  if (!all(names(weights) %in% colnames(m))) stop("weight names do not match the component panel")
  idx <- as.numeric(m[, names(weights), drop = FALSE] %*% weights)
  names(idx) <- rownames(m)
  idx
}

## y-axis breaks: pick the candidate step
## closest to a six-break grid and round the limits to it (DescTools::RoundTo = round(x / gap) * gap).
axis_breaks <- function(x) {
  init <- seq(from = x[1], to = x[2], length.out = 6)
  step <- abs(init[2] - init[1])
  candidates <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10)
  gap <- candidates[which.min(abs(step - candidates))]
  xr <- round(x / gap) * gap
  seq(from = xr[1], to = xr[2], by = gap)
}

## NBER recessions (data/nber_recessions.csv) for shading from `from` on; a recession that
## straddles the first plotted month is clipped to it (as in the published figures).
recession_frame <- function(data_dir, from) {
  nber <- utils::read.csv(file.path(data_dir, "nber_recessions.csv"), stringsAsFactors = FALSE)
  nber$start <- as.Date(nber$start); nber$end <- as.Date(nber$end)
  from <- as.Date(from)
  straddle <- nber$start <= from & nber$end >= from
  nber$start[straddle] <- from
  nber[nber$start >= from, , drop = FALSE]
}

## ---- assembled-index time series (fitted_aci.png, fitted_noaa.png) ----------------------------

## series: named list of numeric vectors named by date (YYYY-MM-DD), in legend order; the
## last element is the headline ACI. dates: the plotting window (sets the x-axis range); a
## series without a value at a date is NA there. Returns the long plot data frame.
series_plot_frame <- function(series, dates) {
  plot_df <- do.call(rbind, lapply(names(series), function(nm) {
    data.frame(date = as.Date(dates), model = nm, value = as.numeric(series[[nm]][dates]), stringsAsFactors = FALSE)
  }))
  plot_df$model <- factor(plot_df$model, levels = names(series))
  plot_df$linewidth <- ifelse(plot_df$model == names(series)[length(series)], 0.8, 1)   # headline thinner
  plot_df
}

## plot_df: from series_plot_frame(). colours: named by model (legend order). labels: expression
## vector in the same order. recessions: from recession_frame(), or NULL for no shading (the
## ACI figure in the paper has none).
draw_series_plot <- function(plot_df, colours, labels, recessions = NULL) {
  linetypes <- stats::setNames(rep("solid", length(colours)), names(colours))
  format_2d <- function(x) sprintf("%.2f", x)
  shading <- if (!is.null(recessions))
    geom_rect(data = recessions, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf), fill = "blue2", alpha = 0.05)
  ggplot(data = plot_df, mapping = aes(x = date)) +
    geom_line(mapping = aes(x = date, y = value, color = model, linetype = model),
              linewidth = plot_df$linewidth, na.rm = TRUE) +
    theme_bw() +
    xlab("") + ylab("") +
    scale_x_date(breaks = scales::pretty_breaks(n = 8), date_labels = "%Y", expand = c(0.01, 0.01)) +
    geom_hline(yintercept = 0, linewidth = 0.6, colour = "grey29") +
    shading +
    scale_color_manual(name = "", breaks = names(colours), values = colours, labels = labels) +
    scale_linetype_manual(name = "", breaks = names(linetypes), values = linetypes, labels = labels, guide = "none") +
    scale_y_continuous(breaks = function(x) axis_breaks(x), labels = format_2d) +
    theme_publication() +
    theme(legend.text = element_text(size = 25),
          legend.key.size = grid::unit(3, "line"),
          legend.margin = margin(5, 5, 5, 5),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.position = "inside",
          legend.position.inside = c(0.20, 0.87),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.key = element_rect(fill = "transparent", colour = NA),
          legend.spacing.x = grid::unit(0.3, "cm"),
          axis.text.x = element_text(size = 25, angle = 0),
          axis.text.y = element_text(size = 25),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(face = "bold", colour = "white", size = 30),
          strip.background = element_rect(colour = "black", fill = "black"),
          plot.margin = grid::unit(c(2, 7, 2, -9), "mm"),
          aspect.ratio = 0.32) +
    guides(col = guide_legend(title.position = "top", nrow = 2, byrow = TRUE,
                              override.aes = list(linewidth = 1.5, linetype = "solid")))
}

## ---- component weights, ACI national + regional (weights_comp_aci.png) ------------------------

## w_regional: 34 regional weights (names = component names); w_national: 5 national weights
## (names = variable names); labels_regional: data/labels_aci_regional.csv.
draw_weights_comp_aci <- function(w_regional, w_national, labels_regional) {
  plot_df <- data.frame(value = unname(w_regional[labels_regional$name]) * 100,
                        names = labels_regional$name, region = labels_regional$region,
                        variable = labels_regional$variable, stringsAsFactors = FALSE)
  ## the WP90 block comes first, as in the paper (rows 28:34, then 1:27)
  variable_order <- c("WP90", "Sea.Level", "Rx5Day", "T10", "T90")
  plot_df <- plot_df[order(match(plot_df$variable, variable_order)), ]
  plot_df$variable <- factor(plot_df$variable, levels = unique(plot_df$variable))
  plot_df$id <- seq_len(nrow(plot_df))

  ## one block per variable: bar range, label position, national weight, equal weight
  var_groups <- data.frame(variable = levels(plot_df$variable), stringsAsFactors = FALSE)
  var_groups$start <- as.numeric(tapply(plot_df$id, plot_df$variable, min)[var_groups$variable])
  var_groups$end   <- as.numeric(tapply(plot_df$id, plot_df$variable, max)[var_groups$variable])
  var_groups$mid   <- (var_groups$start + var_groups$end) / 2
  var_groups$agg_value <- unname(w_national[var_groups$variable]) * 100
  var_groups$Official  <- 100 / nrow(var_groups)

  all_fill_colors <- c(ALA = "#D44233", CWP = "#197A81", SEA = "#782F26", SWP = "#243835",
                       CEA = "#3C68BB", MID = "#EFAE0E", SPL = "#FF8500")

  ## the national weights are drawn on a secondary axis, rescaled to the regional maximum
  max_left <- max(plot_df$value, var_groups$Official, na.rm = TRUE)
  y_top <- max_left * 1.02
  scale_factor <- max(plot_df$value, var_groups$Official, na.rm = TRUE) / max(var_groups$agg_value, na.rm = TRUE)
  var_groups$agg_value_scaled <- var_groups$agg_value * scale_factor

  gg_coef <- ggplot() +
    geom_rect(data = var_groups,
              aes(xmin = start - 0.5, xmax = end + 0.5, ymin = 0, ymax = agg_value_scaled, linetype = "National"),
              fill = "#001A43", color = "#001A43", linewidth = 0.0, alpha = 0.10) +
    geom_bar(data = plot_df[which(!is.na(plot_df$value)), ],
             mapping = aes(x = id, y = value, fill = region),
             stat = "identity", alpha = 0.85, width = 0.8) +
    theme_bw() + xlab("") + ylab("") +
    scale_y_continuous(name = "Weights (in %)", expand = c(0, 0), limits = c(0, y_top),
                       sec.axis = sec_axis(~ . / scale_factor, name = NULL)) +
    ## legend key of the national weights, drawn as a rectangle at fixed data coordinates
    annotate("rect", xmin = 14, xmax = 14.8, ymin = 0.82 * y_top, ymax = 0.91 * y_top,
             fill = "#001A43", color = NA, alpha = 0.10) +
    annotate("text", x = Inf, y = Inf,
             label = "bold(MRCI[scriptstyle(comps)]^{scriptstyle(Nation)}) ~ 'weights (in %) \u2192'",
             parse = TRUE, hjust = 1.05, vjust = 1.5, size = 11.5, color = "black") +
    scale_x_continuous(breaks = var_groups$mid, labels = var_groups$variable, expand = c(0.01, 0.01)) +
    theme(strip.text = element_text(face = "bold", colour = "white", size = 35),
          strip.background = element_rect(colour = "black", fill = "black"),
          plot.margin = grid::unit(c(4, 9, -20, 9), "mm"),
          aspect.ratio = 0.54,
          axis.text.x = element_text(size = 32.6, angle = 0, hjust = 0.5, vjust = -0.5),
          axis.ticks.x = element_line(),
          axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          axis.text.y = element_text(size = 32.6),
          axis.title.y = element_text(size = 32.6, angle = 90, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.8)) +
    scale_linetype_manual(name = "", values = c(National = "solid"),
                          labels = expression(MRCI[scriptstyle(comps)]^{scriptstyle(Nation)})) +
    scale_fill_manual(name = expression(MRCI[scriptstyle(comps)]^{scriptstyle(Region)} * " :"),
                      values = all_fill_colors) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 32.6),
          legend.title = element_text(size = 32.6),
          legend.key.size = grid::unit(1.5, "line"),
          legend.spacing.x = grid::unit(1.5, "cm"),
          legend.margin = margin(12, 0, 0, 0),
          legend.box.margin = margin(5, 0, 0, 0)) +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black")
  gg_coef
}

## ---- component weights, NOAA aggregated to region x variable (weights_comp_noaa.png) ---------

## w: 82 NOAA weights in panel order; labels_noaa: data/labels_noaa.csv (region_agg and
## variable_agg give the aggregation of Tables C.1-C.2).
draw_weights_comp_noaa <- function(w, labels_noaa) {
  cells <- data.frame(region = labels_noaa$region_agg, variable = labels_noaa$variable_agg,
                      weight = as.numeric(w)[labels_noaa$index], stringsAsFactors = FALSE)
  agg <- stats::aggregate(weight ~ region + variable, data = cells, FUN = sum)

  variable_order <- c("Heat.DD", "Palmer.Z", "Precip.", "Temp.", "Others")
  region_order   <- c("Up.MidW.", "Ohio.V.", "N.West", "N.Rockies", "N.East", "Rest")
  plot_df <- data.frame(value = agg$weight * 100,
                        variable = factor(agg$variable, levels = variable_order),
                        region = factor(agg$region, levels = region_order))
  plot_df <- plot_df[order(plot_df$variable, plot_df$region), ]
  plot_df$id <- seq_len(nrow(plot_df))
  var_groups <- data.frame(variable = levels(plot_df$variable), stringsAsFactors = FALSE)
  var_groups$mid <- as.numeric(tapply(plot_df$id, plot_df$variable, function(i) mean(range(i)))[var_groups$variable])

  region_colors <- c("Up.MidW." = "#D44233", "Ohio.V." = "#197A81", "N.West" = "#FF8500",
                     "N.Rockies" = "#3C68BB", "N.East" = "#EFAE0E", "Rest" = "#782F26")
  minor <- c(6.5, 12.5, 18.5, 24.5)          # dashed separators between the variable blocks

  ggplot() +
    geom_bar(data = plot_df, aes(x = id, y = value, fill = region),
             stat = "identity", alpha = 0.85, width = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    scale_y_continuous(name = "Weights (in %)", expand = c(0, 0), limits = c(0, max(plot_df$value) * 1.02)) +
    scale_x_continuous(breaks = var_groups$mid, labels = var_groups$variable,
                       expand = c(0.01, 0.01), minor_breaks = minor) +
    scale_fill_manual(name = expression(MRCI[scriptstyle(comps)]^{scriptstyle(NOAA)} * " :"),
                      values = region_colors) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(color = "gray45", linewidth = 1, linetype = 2),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 31.2, angle = 90, face = "bold"),
          axis.text.x = element_text(size = 31.2, angle = 0, hjust = 0.5, vjust = 0),
          axis.text.y = element_text(size = 31.2),
          axis.ticks.x = element_line(),
          plot.margin = grid::unit(c(4, 9, -20, 9), "mm"),
          aspect.ratio = 0.54,
          legend.position = "bottom",
          legend.text = element_text(size = 31.2),
          legend.title = element_text(size = 31.2),
          legend.key.size = grid::unit(1.2, "line"),
          legend.key.spacing.y = grid::unit(5, "pt"),
          legend.spacing.x = grid::unit(3, "cm"),
          legend.margin = margin(3, 0, 0, 0)) +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE))
}

## ---- rank weights (weights_rank_aci.png, weights_rank_noaa.png) -------------------------------

## Displayed rank profiles: |w| normalised, then smoothed.
## ACI regional: trapezoidal kernel over eight ranks (left-aligned 4-rank then 5-rank means).
smooth_rank_weights_trapezoid <- function(w) {
  w <- abs(as.numeric(w)); w <- w / sum(w, na.rm = TRUE)
  w <- zoo::rollmean(w, k = 4, align = "left")
  w <- zoo::rollmean(w, k = 5, align = "left")
  w / sum(w)
}
## NOAA: sliding average over two adjacent ranks.
smooth_rank_weights_pair <- function(w) {
  w <- abs(as.numeric(w)); w <- w / sum(w, na.rm = TRUE)
  w <- (w[-length(w)] + w[-1]) / 2
  w / sum(w)
}

## w_smoothed: from one of the smoothers. legend_label: expression. text_size / aspect_ratio:
## 36.4 / 0.572 (ACI) and 35.6 / 0.575 (NOAA), as in the published figures.
draw_weights_rank <- function(w_smoothed, legend_label, text_size = 36.4, aspect_ratio = 0.572) {
  plot_df <- data.frame(value = w_smoothed * 100, rank = seq_along(w_smoothed), group = "rank_weights",
                        stringsAsFactors = FALSE)
  quintiles <- round(stats::quantile(plot_df$rank, c(0.2, 0.4, 0.6, 0.8, 1)))   # x labels 20% .. 100%
  right_margin <- 8.6 + 6

  ggplot(data = plot_df) +
    geom_area(aes(x = rank, y = value, fill = group, group = group), alpha = 0.55, position = "identity") +
    geom_line(aes(x = rank, y = value, colour = group, group = group), linewidth = 1, show.legend = FALSE) +
    geom_point(aes(x = rank, y = value, colour = group, group = group), size = 1, show.legend = FALSE) +
    theme_bw() + xlab("") + ylab("") +
    scale_y_continuous(name = NULL, expand = c(0, 0),
                       limits = c(0, max(plot_df$value) + (max(plot_df$value) / 100)),
                       breaks = scales::pretty_breaks(n = 5),
                       sec.axis = dup_axis(name = "Weights (in %)")) +          # right-hand label only
    theme_publication() +
    scale_x_continuous(breaks = quintiles, labels = names(quintiles), expand = c(0, 0)) +
    scale_fill_manual(name = "", values = c(rank_weights = "#001330"), labels = legend_label,
                      guide = guide_legend(direction = "horizontal", title.position = "top", title.hjust = 0.5)) +
    scale_colour_manual(name = "", values = c(rank_weights = "#001330"), guide = "none") +
    theme(strip.text = element_text(face = "bold", colour = "white", size = 40),
          strip.background = element_rect(colour = "black", fill = "black"),
          plot.margin = grid::unit(c(4, right_margin, 5, 5), "mm"),
          aspect.ratio = aspect_ratio,
          axis.title.y.left = element_blank(),
          axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank(),
          axis.title.y.right = element_text(size = text_size, face = "bold"),
          axis.text.y.right = element_text(size = text_size),
          axis.ticks.y.right = element_line(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = text_size, angle = 0, vjust = -0.5),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "bottom",
          legend.text = element_text(size = text_size),
          legend.key.size = grid::unit(2, "line"),
          legend.spacing.x = grid::unit(12, "cm"),
          legend.margin = margin(-5, 0, 0, 30, unit = "mm"),
          legend.title = element_text(size = 10, face = "bold")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 0.55)))
}
