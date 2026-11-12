## =============================================================================
## code/funcs/load_data.R -- read the shipped cleaned panels and build the estimation
## sample for one specification.
##
## The shipped files in data/ are seasonally adjusted LEVEL panels (from 1961)
## plus, for auditability, the TRANSFORMED panels exactly as the estimator uses
## them. This loader re-derives the transformed panels from the level panels and
## stops if they do not match the shipped transformed files (see the data
## dictionary, data/README.md, for every formula).
##
## Returns a list with the trimmed estimation sample (covid months removed),
## the untrimmed calendar and gap dummies, and the pre-trim panels.
## =============================================================================

## ---- helpers ------------------------------------------------------------------------------

read_panel <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- x[[1]]; x <- x[, -1, drop = FALSE]
  x <- as.data.frame(lapply(x, as.numeric), check.names = FALSE, row.names = rn)
  x
}

## month-on-month log difference in percent; first row is dropped
log_diff_pct <- function(x) {
  n <- nrow(x)
  out <- (log(x[2:n, , drop = FALSE]) - log(x[1:(n - 1), , drop = FALSE])) * 100
  rownames(out) <- rownames(x)[2:n]
  out
}

## Standardise each column with the monthly mean and standard deviation over the
## reference years (the Actuaries Climate Index convention, 1961-1990).
standardize_monthly <- function(x, ref_years = 1961:1990) {
  d  <- as.Date(rownames(x))
  yr <- as.integer(format(d, "%Y")); mo <- as.integer(format(d, "%m"))
  out <- x
  for (m in 1:12) {
    ref <- which(mo == m & yr %in% ref_years); rows <- which(mo == m)
    mu <- vapply(x[ref, , drop = FALSE], mean, numeric(1), na.rm = TRUE)
    sg <- vapply(x[ref, , drop = FALSE], stats::sd, numeric(1), na.rm = TRUE)
    out[rows, ] <- sweep(sweep(as.matrix(x[rows, , drop = FALSE]), 2, mu, "-"), 2, sg, "/")
  }
  out
}

## Row-wise sort of a panel (rank space): value in column j is the j-th order statistic.
rank_sort <- function(x) {
  m <- as.matrix(x)
  sorted <- t(apply(m, 1, sort, na.last = FALSE))
  perm   <- t(apply(m, 1, order, na.last = FALSE))
  dimnames(sorted) <- list(rownames(m), paste0("rank_", seq_len(ncol(m))))
  dimnames(perm)   <- dimnames(sorted)
  list(sorted = as.data.frame(sorted, check.names = FALSE), perm = perm)
}

## NOAA per-series stationarity transformation (rule chosen by a KPSS test when the
## panel was built; the rules are shipped in data/noaa_transformation_rules.csv).
transform_noaa_series <- function(x, rule) {
  n <- length(x); out <- rep(NA_real_, n)
  if (rule == "demean") {
    out <- x - mean(stats::na.omit(x)); out[1] <- NA_real_
  } else if (rule == "logdiff") {
    x2 <- x; x2[x2 <= 0] <- NA_real_
    out[2:n] <- log(x2[2:n]) - log(x2[1:(n - 1)])
  } else if (rule == "diff") {
    out[2:n] <- x[2:n] - x[1:(n - 1)]
  } else if (rule == "zero") {
    out <- rep(0, n)
  } else stop("unknown NOAA transformation rule: ", rule)
  out
}

assert_same <- function(a, b, what, tol = 1e-10) {
  a <- as.matrix(a); b <- as.matrix(b)
  ok <- isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE)) &&
        identical(dim(a), dim(b)) && identical(rownames(a), rownames(b))
  if (!ok) stop("re-derived ", what, " does not match the shipped transformed file")
  invisible(TRUE)
}

## ---- climate panels -----------------------------------------------------------------------

## National ACI: 6 seasonally adjusted components -> MoM log-diff (%) -> monthly standardisation
## on 1961-1990 -> T10 sign flipped -> CDD dropped. The headline (equal-weight) ACI is the
## row mean of the same five standardised components.
build_aci_national <- function(data_dir, check = TRUE) {
  lvl <- read_panel(file.path(data_dir, "aci_national_level.csv"))
  std <- standardize_monthly(log_diff_pct(lvl))
  std$T10 <- -std$T10
  comps <- std[, setdiff(colnames(std), "CDD"), drop = FALSE]
  headline <- rowMeans(comps); names(headline) <- rownames(comps)
  if (check) {
    assert_same(comps, read_panel(file.path(data_dir, "aci_national_transformed.csv")), "national ACI panel")
    assert_same(matrix(headline, ncol = 1, dimnames = list(names(headline), "ACI")),
                read_panel(file.path(data_dir, "aci_headline.csv")), "headline ACI")
  }
  list(comps = comps, headline = headline)
}

## Regional ACI: 34 seasonally adjusted region-components -> absolute value -> MoM log-diff (%)
## -> monthly standardisation -> all T10.* columns sign flipped.
build_aci_regional <- function(data_dir, check = TRUE) {
  lvl <- read_panel(file.path(data_dir, "aci_regional_level.csv"))
  std <- standardize_monthly(log_diff_pct(abs(lvl)))
  t10 <- grep("T10", colnames(std), value = TRUE); std[, t10] <- -std[, t10]
  if (check) assert_same(std, read_panel(file.path(data_dir, "aci_regional_transformed.csv")), "regional ACI panel")
  std
}

## NOAA panel: 82 seasonally adjusted series -> per-series demean / log-diff / diff -> first
## month dropped -> internal gaps linearly interpolated -> monthly standardisation -> non-finite
## values set to zero.
build_noaa <- function(data_dir, check = TRUE) {
  lvl   <- read_panel(file.path(data_dir, "noaa_level.csv"))
  rules <- utils::read.csv(file.path(data_dir, "noaa_transformation_rules.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(rules$name, colnames(lvl)))
  st <- mapply(function(cn, rule) transform_noaa_series(lvl[[cn]], rule), rules$name, rules$rule, SIMPLIFY = TRUE)
  st <- st[-1, , drop = FALSE]; rownames(st) <- rownames(lvl)[-1]
  st[!is.finite(st)] <- NA
  stopifnot(all(colMeans(is.na(st)) <= 0.075), all(apply(st, 2, stats::sd, na.rm = TRUE) > 0))
  st <- apply(st, 2, function(v) zoo::na.approx(v, na.rm = FALSE)); rownames(st) <- rownames(lvl)[-1]
  stopifnot(!anyNA(st))
  z <- as.matrix(standardize_monthly(as.data.frame(st, check.names = FALSE)))
  z[!is.finite(z)] <- 0
  z <- as.data.frame(z, check.names = FALSE)
  if (check) assert_same(z, read_panel(file.path(data_dir, "noaa_transformed.csv")), "NOAA panel")
  z
}

## ---- macro block ----------------------------------------------------------------------------

## FRED-MD levels -> 100*dlog (INDPRO, PCEPI, EXUSUKx, HOUST), first difference (UNRATE),
## level (SHADOWRATE). The first observation keeps its level value; it is discarded by the
## 1975 sample start.
build_macro <- function(data_dir, log_diff = c("INDPRO", "PCEPI", "EXUSUKx", "HOUST"), diff = "UNRATE") {
  m <- read_panel(file.path(data_dir, "macro_fredmd_2025-12.csv"))
  n <- nrow(m)
  m[2:n, log_diff] <- (log(m[2:n, log_diff]) - log(m[1:(n - 1), log_diff])) * 100
  m[2:n, diff]     <- m[2:n, diff] - m[1:(n - 1), diff]
  m
}

## ---- assemble one specification ---------------------------------------------------------------

## spec: list(level = 1|2|3, rank = TRUE|FALSE). cfg: list with var_order (climate first), sample
## start/end, covid dates, end_test_date.
load_specification <- function(spec, data_dir, cfg, check = TRUE) {
  level <- spec$level; rank <- isTRUE(spec$rank)
  clim_name <- cfg$climate_name                       # column name of the climate variable
  var_order <- c(clim_name, cfg$macro_order)
  macro_names <- cfg$macro_order

  nat <- build_aci_national(data_dir, check)
  comps <- switch(as.character(level), "1" = nat$comps, "2" = build_aci_regional(data_dir, check), "3" = build_noaa(data_dir, check))
  perm <- NULL
  if (rank) { rs <- rank_sort(comps); comps <- rs$sorted; perm <- rs$perm }

  macro <- build_macro(data_dir)
  macro[[clim_name]] <- unname(nat$headline[rownames(macro)])   # equal-weight benchmark for all levels
  macro <- macro[, var_order]

  ## common sample and trimming used by the saved estimates
  common <- intersect(rownames(macro), rownames(comps))
  end_date   <- min(cfg$end_date, max(common)); start_date <- max(cfg$start_date, min(common))
  end_test   <- cfg$end_test_date
  if (end_test < start_date || end_test > end_date) end_test <- end_date
  common <- common[which(common == start_date):which(common == end_date)]
  common_full <- common
  covid_in <- cfg$covid_dates[cfg$covid_dates %in% common_full]
  gap_dummy <- NULL
  if (length(covid_in)) {
    gap_dummy <- matrix(0, length(common_full), length(covid_in),
                        dimnames = list(common_full, paste0("d_", gsub("-", "", covid_in))))
    for (j in seq_along(covid_in)) gap_dummy[match(covid_in[j], common_full), j] <- 1
    common <- common[!common %in% cfg$covid_dates]
  }
  raw_bench <- macro[rownames(macro) %in% common_full, ]
  raw_comps <- comps[rownames(comps) %in% common_full, , drop = FALSE]
  bench <- macro[rownames(macro) %in% common, ]                       # all variables incl. climate benchmark
  mac   <- macro[rownames(macro) %in% common, macro_names]            # macro variables only
  cmp   <- comps[rownames(comps) %in% common, , drop = FALSE]
  if (!is.null(perm)) perm <- perm[rownames(perm) %in% common, , drop = FALSE]

  list(spec = spec, var_order = var_order, climate_name = clim_name,
       dates = common, dates_full = common_full, gap_dummy = gap_dummy,
       start_date = start_date, end_date = end_date, end_test_date = end_test,
       macro = mac, bench = bench, comps = cmp, comps_perm = perm,
       raw_bench = raw_bench, raw_comps = raw_comps,
       shrink_to = rep(1 / ncol(cmp), ncol(cmp)),
       pre = list(macro = macro, comps = comps, headline = nat$headline))
}
