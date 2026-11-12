## =============================================================================
## Numerical helpers for Appendix A/B. They draw no figures; the figure and
## summary scripts use them. All are pure except the
## assembled-index helpers at the end of this file, which read the shipped
## seasonally adjusted level panels from data/ (see smoothed_level_panel).
## =============================================================================

usable_boot_weights <- function(task1) {
  W <- as.matrix(task1$boot$weights)
  keep <- as.logical(task1$boot$converged)
  keep[is.na(keep)] <- FALSE
  keep <- keep & stats::complete.cases(W)
  W[keep, , drop = FALSE]
}

appendix_groups <- function(id, point, labels_regional, labels_noaa) {
  if (id == "comp_national") return(names(point))
  if (id == "comp_regional") return(labels_regional$variable)
  if (id == "comp_noaa") return(labels_noaa$variable_agg)
  stop("component-family groups are not defined for ", id)
}

aggregate_weight_rows <- function(W, groups, levels = unique(groups)) {
  W <- as.matrix(W)
  out <- vapply(levels, function(g) rowSums(abs(W[, groups == g, drop = FALSE])), numeric(nrow(W)))
  if (is.null(dim(out)))
    out <- if (nrow(W) == 1L) matrix(out, nrow = 1L) else matrix(out, ncol = 1L)
  colnames(out) <- levels
  out
}

smooth_rank_profile <- function(w, id) {
  w <- abs(as.numeric(w)); w <- w / sum(w)
  if (id == "rank_regional") {
    w <- zoo::rollmean(w, 4L, align = "left")
    w <- zoo::rollmean(w, 5L, align = "left")
  } else if (id == "rank_noaa") {
    w <- (w[-length(w)] + w[-1L]) / 2
  }
  w / sum(w)
}

bootstrap_weight_statistics <- function(id, point_weights, task1, labels_regional, labels_noaa) {
  W <- usable_boot_weights(task1)
  probs <- c(.05, .16, .50, .84, .95)
  if (startsWith(id, "comp_")) {
    groups <- appendix_groups(id, point_weights, labels_regional, labels_noaa)
    lev <- if (id == "comp_national") names(point_weights) else unique(groups)
    point <- aggregate_weight_rows(matrix(point_weights, nrow = 1L), groups, lev)[1L, ]
    draws <- aggregate_weight_rows(W, groups, lev)
  } else {
    point <- smooth_rank_profile(point_weights, id)
    draws <- t(apply(W, 1L, smooth_rank_profile, id = id))
    lev <- paste0("rank_", seq_along(point))
    colnames(draws) <- lev; names(point) <- lev
  }
  q <- t(apply(draws, 2L, stats::quantile, probs = probs, names = FALSE, na.rm = TRUE))
  data.frame(id = id, group = lev, point = as.numeric(point),
             q05 = q[, 1L], q16 = q[, 2L], q50 = q[, 3L], q84 = q[, 4L], q95 = q[, 5L],
             n_draws = nrow(draws), stringsAsFactors = FALSE)
}

trailing_mean_matrix <- function(x, k) {
  x <- as.matrix(x)
  if (k <= 1L) return(x)
  out <- apply(x, 2L, function(z) zoo::rollmean(z, k, fill = NA_real_, align = "right"))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  rownames(out) <- rownames(x)
  out
}

bootstrap_index_statistics <- function(id, raw_comps, point_weights, task1) {
  W <- usable_boot_weights(task1)
  X <- as.matrix(raw_comps)
  point <- X %*% as.numeric(point_weights)
  paths <- X %*% t(W)
  k <- if (id == "rank_noaa") 12L else if (startsWith(id, "rank_")) 24L else 1L
  point <- trailing_mean_matrix(point, k)[, 1L]
  paths <- trailing_mean_matrix(paths, k)
  q <- t(apply(paths, 1L, stats::quantile, probs = c(.05, .16, .50, .84, .95),
               names = FALSE, na.rm = TRUE))
  ok <- is.finite(point) & is.finite(q[, 2L]) & is.finite(q[, 4L])
  data.frame(id = id, date = as.Date(rownames(X)), point = as.numeric(point),
             q05 = q[, 1L], q16 = q[, 2L], q50 = q[, 3L], q84 = q[, 4L], q95 = q[, 5L],
             covered68 = ifelse(ok, point >= q[, 2L] & point <= q[, 4L], NA),
             covered90 = ifelse(ok, point >= q[, 1L] & point <= q[, 5L], NA),
             stringsAsFactors = FALSE)
}

expanding_weight_calendar <- function(W_mat, smooth = TRUE,
                                      start = "2005-10-01", end = "2025-03-01",
                                      covid_start = "2020-01-01", covid_end = "2021-06-01") {
  W_mat <- as.matrix(W_mat)
  dates <- seq.Date(as.Date(start), as.Date(end), by = "month")
  covid <- seq.Date(as.Date(covid_start), as.Date(covid_end), by = "month")
  keep <- !dates %in% covid
  if (sum(keep) != nrow(W_mat))
    stop("expanding matrix has ", nrow(W_mat), " rows; expected ", sum(keep))
  out <- matrix(NA_real_, length(dates), ncol(W_mat),
                dimnames = list(as.character(dates), colnames(W_mat)))
  out[keep, ] <- W_mat
  anchor <- match(as.Date("2019-12-01"), dates)
  out[!keep, ] <- matrix(out[anchor, ], sum(!keep), ncol(out), byrow = TRUE)
  if (smooth) {
    out <- apply(out, 2L, function(z) zoo::rollapply(z, 12L, mean, partial = TRUE,
                                                     align = "center", na.rm = TRUE))
    if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
    rownames(out) <- as.character(dates)
  }
  out
}

expanding_family_statistics <- function(id, W_mat, point_weights, labels_regional, labels_noaa) {
  cal <- expanding_weight_calendar(W_mat, smooth = TRUE)
  groups <- appendix_groups(id, point_weights, labels_regional, labels_noaa)
  lev <- if (id == "comp_national") names(point_weights) else unique(groups)
  shares <- aggregate_weight_rows(cal, groups, lev)
  shares <- shares / rowSums(shares)
  data.frame(id = id, group = lev, start = shares[1L, ], end = shares[nrow(shares), ],
             mean = colMeans(shares), min = apply(shares, 2L, min),
             min_date = as.Date(rownames(shares)[apply(shares, 2L, which.min)]),
             stringsAsFactors = FALSE)
}

rank_usage_noaa <- function(rank_weights, component_weights, permutation, labels_noaa) {
  K <- length(rank_weights); permutation <- as.matrix(permutation)
  occ <- matrix(0, K, K)
  for (j in seq_len(K)) occ[, j] <- tabulate(permutation[, j], nbins = K) / nrow(permutation)
  usage <- K * as.numeric(occ %*% as.numeric(rank_weights)) - 1
  cells <- stats::aggregate(usage,
    list(region = labels_noaa$region_agg, family = labels_noaa$variable_agg), mean)
  names(cells)[3L] <- "usage"
  fam_usage <- stats::aggregate(cells$usage, list(family = cells$family), mean)
  names(fam_usage)[2L] <- "usage"
  comp_cells <- stats::aggregate(as.numeric(component_weights),
    list(region = labels_noaa$region_agg, family = labels_noaa$variable_agg), sum)
  names(comp_cells)[3L] <- "weight"
  comp_family <- stats::aggregate(comp_cells$weight, list(family = comp_cells$family), sum)
  names(comp_family)[2L] <- "weight"
  ord <- match(fam_usage$family, comp_family$family)
  list(component_usage = usage, cells = cells, family = fam_usage,
       component_family = comp_family,
       spearman_family = stats::cor(fam_usage$usage, comp_family$weight[ord], method = "spearman"))
}

## ---- Appendix A.2: the plotted assembled-index envelope ---------------------------------------
## The published index panels are drawn on FIVE-YEAR moving averages of the seasonally adjusted
## level panels, not on the monthly transformed panel that the estimator uses.  The construction
## below is the one used for the published figures.  Each level series is interpolated over
## internal gaps and averaged over a 60-month window labelled by the LAST month of that window;
## T10 columns are sign-flipped so that every component increases with adverse climate stress.
##
## The two panels differ in the order of standardisation: the ACI panels follow the convention of
## the published index, which standardises the MONTHLY series on the 1961-1990 moments and smooths
## afterwards, so the smoothed series keeps the smaller variance of a five-year average.  The NOAA
## panel standardises AFTER smoothing, which rescales the smoothed series to unit variance.
##
## These helpers need read_panel(), standardize_monthly() and the build_*() panel constructors
## from code/funcs/load_data.R, which the calling script sources.

smoothed_level_panel <- function(panel, data_dir, k = 60L) {
  file <- switch(panel,
                 national = "aci_national_level.csv",
                 regional = "aci_regional_level.csv",
                 noaa     = "noaa_level.csv",
                 stop("unknown panel: ", panel))
  standardize_first <- panel %in% c("national", "regional")
  lvl <- read_panel(file.path(data_dir, file))
  X <- apply(as.matrix(lvl), 2L, function(z) zoo::na.approx(z, na.rm = FALSE))
  rownames(X) <- rownames(lvl)
  if (standardize_first)
    X <- as.matrix(standardize_monthly(as.data.frame(X, check.names = FALSE)))
  sm <- zoo::rollmean(X, k)
  rownames(sm) <- rownames(lvl)[k:nrow(lvl)]
  sm <- as.data.frame(sm, check.names = FALSE)
  if (!standardize_first) sm <- standardize_monthly(sm)
  t10 <- grep("T10", colnames(sm), value = TRUE)
  if (length(t10)) sm[, t10] <- -sm[, t10]
  sm
}

## Equal-weight headline ACI on the same smoothed, standardised basis: the row mean of all six
## national components (consecutive dry days included, T10 already flipped).  It is the green
## reference line in every panel, for the NOAA specifications as well.
smoothed_aci_headline <- function(data_dir, k = 60L) {
  sm <- smoothed_level_panel("national", data_dir, k)
  stats::setNames(rowMeans(sm, na.rm = TRUE), rownames(sm))
}

## Constituent panel in component space: as above, without consecutive dry days.
smoothed_component_panel <- function(panel, data_dir, k = 60L) {
  sm <- smoothed_level_panel(panel, data_dir, k)
  as.matrix(sm[, setdiff(colnames(sm), "CDD"), drop = FALSE])
}

## Constituent panel in rank space.  For the regional and NOAA panels the smoothed columns are
## reordered, date by date, by the rank positions of the UNSMOOTHED transformed panel (the object
## the rank weights were estimated on); the national panel is sorted row-wise instead.
smoothed_rank_panel <- function(panel, data_dir, k = 60L) {
  sm <- smoothed_component_panel(panel, data_dir, k)
  if (panel == "national") {
    out <- t(apply(sm, 1L, sort, na.last = FALSE))
    rownames(out) <- rownames(sm)
  } else {
    tr <- switch(panel,
                 regional = build_aci_regional(data_dir),
                 noaa     = build_noaa(data_dir))
    perm <- t(apply(as.matrix(tr), 1L, order, na.last = FALSE))
    rownames(perm) <- rownames(tr)
    keep <- intersect(rownames(sm), rownames(perm))
    sm <- sm[keep, , drop = FALSE]; perm <- perm[keep, , drop = FALSE]
    out <- t(vapply(seq_along(keep), function(i) as.numeric(sm[i, perm[i, ]]), numeric(ncol(sm))))
    rownames(out) <- keep
  }
  colnames(out) <- paste0("rank_", seq_len(ncol(out)))
  out
}

## Point index and the envelope of the converged bootstrap weight vectors, both pushed through the
## construction above.  Rank-space paths carry the extra right-aligned smoothing of the paper
## (24 months for the ACI panels, 12 for NOAA); component-space paths need none, the 60-month
## average of the input panel already provides it.  The COVID months stay in the plotting window.
bootstrap_index_paths <- function(id, spec, data_dir, point_weights, task1,
                                  start = "1975-01-01", end = "2025-05-01") {
  panel <- c("national", "regional", "noaa")[spec$level]
  X <- if (isTRUE(spec$rank)) smoothed_rank_panel(panel, data_dir)
       else smoothed_component_panel(panel, data_dir)
  W <- usable_boot_weights(task1)
  if (ncol(X) != length(point_weights) || ncol(X) != ncol(W))
    stop("smoothed panel has ", ncol(X), " columns; weights have ",
         length(point_weights), " / ", ncol(W))
  point <- X %*% as.numeric(point_weights)
  paths <- X %*% t(W)
  k <- if (!isTRUE(spec$rank)) 1L else if (spec$level == 3L) 12L else 24L
  point <- trailing_mean_matrix(point, k)[, 1L]
  paths <- trailing_mean_matrix(paths, k)
  q <- t(apply(paths, 1L, stats::quantile, probs = c(.05, .16, .50, .84, .95),
               names = FALSE, na.rm = TRUE))
  bench <- smoothed_aci_headline(data_dir, k = 60L)
  out <- data.frame(id = id, date = as.Date(rownames(X)), point = as.numeric(point),
                    benchmark = as.numeric(bench[rownames(X)]),
                    q05 = q[, 1L], q16 = q[, 2L], q50 = q[, 3L], q84 = q[, 4L], q95 = q[, 5L],
                    n_draws = nrow(W), stringsAsFactors = FALSE)
  out[out$date >= as.Date(start) & out$date <= as.Date(end), , drop = FALSE]
}

## ---- Summary numbers (code/06_summary_numbers.R) ----------------------------------------------
## Additions used only by the summary script; the figure scripts do not call them.

## Size-neutral rank-implied usage U_m = d * sum_r w_r p_mr - 1 (Appendix B, Equation B.2) for any
## rank-space panel; p_mr is the share of dates on which component m occupies rank r.
rank_usage_components <- function(rank_weights, permutation) {
  K <- length(rank_weights); permutation <- as.matrix(permutation)
  occ <- matrix(0, K, K)
  for (j in seq_len(K)) occ[, j] <- tabulate(permutation[, j], nbins = K) / nrow(permutation)
  K * as.numeric(occ %*% as.numeric(rank_weights)) - 1
}

## Percentile position (j - 0.5) / K of rank j and the share of the rank weights between two
## percentiles. Displayed (smoothed) profiles of length n are labelled at (i - 1) / (n - 1), the
## convention of the x-axis of weights_rank_*.png (quantiles of the plotted index 1..n).
rank_percentile <- function(K) (seq_len(K) - 0.5) / K
plotted_percentile <- function(n) (seq_len(n) - 1) / (n - 1)
rank_mass <- function(w, lo, hi) {
  w <- abs(as.numeric(w)); q <- rank_percentile(length(w))
  sum(w[q >= lo & q <= hi]) / sum(w)
}

## Share of dates on which the point index lies inside the pointwise 68% and 90% bands of
## bootstrap_index_paths() (Appendix A.2, Figure A.2).
index_band_coverage <- function(paths) {
  p <- paths[is.finite(paths$point) & is.finite(paths$q05) & is.finite(paths$q95), , drop = FALSE]
  c(n = nrow(p), cover68 = mean(p$point >= p$q16 & p$point <= p$q84),
    cover90 = mean(p$point >= p$q05 & p$point <= p$q95))
}

## Expanding-window shares as plotted in Figure A.5: COVID months carried forward, centred partial
## 12-month mean, then aggregated to families (component space) or rank quintiles (rank space).
rank_quintile_groups <- function(K) paste0("Q", pmin(5L, ceiling(seq_len(K) * 5 / K)))
expanding_share_matrix <- function(W_mat, groups, levels = unique(groups)) {
  sh <- aggregate_weight_rows(expanding_weight_calendar(W_mat, smooth = TRUE), groups, levels)
  sh / rowSums(sh)
}
