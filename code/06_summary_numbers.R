## =============================================================================
## 06_summary_numbers.R -- Step 6: the numbers quoted in the text of the paper
## (main text and Appendices A-C), computed from the package results and written
## to results/summary_numbers.txt. Each block names the paper location it
## supports. Also cross-checks the NOAA aggregation tables (Appendix C, Tables
## C.1 and C.2) against data/labels_noaa.csv.
##
## Reads  results/point_estimates/ (step 1), results/irf_bands/ (supplied bands),
##        results/bootstrap/ (archived stability results; _recomputed/ with
##        USE_RECOMPUTED=TRUE), data/
## Writes results/summary_numbers.txt
##
## Conventions. Weights are in percent of the index. Impulse-response horizons
## are h = 0..36 with h = 0 the impact month; the paper counts months from 1
## (month = h + 1). Rank j of K sits at percentile (j - 0.5) / K; a displayed
## (smoothed) rank profile of length n is labelled at (i - 1) / (n - 1), as on
## the x-axis of weights_rank_*.png. Point values use the main-text point
## estimates; intervals use the supplied bootstrap draws, as in the figures.
##
## Usage:  Rscript code/06_summary_numbers.R
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "load_data.R"))
source(file.path(ROOT, "code", "funcs", "irf_bands.R"))
source(file.path(ROOT, "code", "funcs", "fig_irf.R"))
source(file.path(ROOT, "code", "funcs", "appendix_summaries.R"))
suppressWarnings(suppressPackageStartupMessages(library(vars)))

out <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); out <<- c(out, s) }
block <- function(title) { if (length(out)) say(""); say("== %s ==", title) }
pct <- function(x, d = 1) {
  z <- 100 * x
  z[abs(z) < 0.5 * 10^(-d)] <- 0
  formatC(z, format = "f", digits = d)
}
shares <- function(x, d = 1) paste(sprintf("%s %s%%", names(x), pct(x, d)), collapse = ", ")
ym <- function(x) substr(as.character(x), 1L, 7L)
## horizons as compact ranges, e.g. 0-3, 7-36
hrange <- function(h) {
  if (!length(h)) return("none")
  s <- split(h, cumsum(c(1, diff(h) != 1)))
  paste(vapply(s, function(x) if (length(x) == 1L) as.character(x) else paste0(x[1], "-", x[length(x)]), ""),
        collapse = ", ")
}

est <- lapply(names(SPECS), function(id)
  readRDS(require_result(file.path(PATHS$estimates, paste0(id, ".rds")), paste0("point estimate for ", id))))
names(est) <- names(SPECS)
bands <- lapply(c("benchmark", names(SPECS)), function(id)
  readRDS(require_result(irf_result_file(id), paste0("IRF band result for ", id))))
names(bands) <- c("benchmark", names(SPECS))
task1 <- lapply(names(SPECS), function(id) {
  t1 <- readRDS(require_result(boot_result_file(id, 1), paste0("task-1 bootstrap result for ", id)))
  names(t1$weights) <- names(est[[id]]$weights)
  t1
})
task2 <- lapply(names(SPECS), function(id)
  readRDS(require_result(boot_result_file(id, 2), paste0("task-2 expanding-window result for ", id))))
names(task1) <- names(task2) <- names(SPECS)
lab_reg  <- read.csv(file.path(PATHS$data, "labels_aci_regional.csv"), stringsAsFactors = FALSE)
lab_noaa <- read.csv(file.path(PATHS$data, "labels_noaa.csv"), stringsAsFactors = FALSE)
comp_ids <- c("comp_national", "comp_regional", "comp_noaa")
rank_ids <- c("rank_national", "rank_regional", "rank_noaa")

## VAR re-fit on the regenerated point estimates (as in step 2), for the sample sizes and to
## confirm that the regenerated weights reproduce the supplied point responses
refit <- lapply(names(bands), function(id) {
  r <- est[[if (id == "benchmark") "comp_national" else id]]
  Y <- cbind(climate = if (id == "benchmark") r$benchmark else r$index, as.matrix(r$macro[, CFG$macro_order]))
  rownames(Y) <- r$dates
  holdout <- length(r$dates) - match(r$end_test_date, r$dates)
  m <- fit_var(Y, r$lags, r$dates, holdout = holdout, use_dummies = id != "benchmark")
  cumulate <- vapply(bands[[id]]$bands, function(b) isTRUE(b$cumulated), logical(1))
  list(irf = irf_matrix(m, "climate", CFG$irf_horizon, cumulate), n_y = nrow(m$y),
       n_resid = nrow(residuals(m)), n_gap = sum(startsWith(colnames(m$datamat), "excl_")))
})
names(refit) <- names(bands)

say("Numbers quoted in the paper, computed from the package results.")
say("Written by code/06_summary_numbers.R. [..] = paper location.")
say("Weights in percent of the index. IRF horizon h = 0..36, h = 0 the impact month; paper month = h + 1.")
say("Appendix inputs: %s pairs-bootstrap and expanding-window results.",
    result_source(boot_result_file("comp_national", 1)))

## ---- data, sample, VAR ----------------------------------------------------------------------
block("[Sec. 1, Sec. 2.4, Sec. 3.1, notes to Figs. 1-2] Panels, sample and VAR")
reg_regions <- unique(lab_reg$region)
say("Constituents K: national %d (%s; CDD excluded), regional %d (%d regions x 5 variables; no Sea.Level in %s), NOAA %d",
    length(est$comp_national$weights), paste(names(est$comp_national$weights), collapse = ", "),
    nrow(lab_reg), length(reg_regions), paste(setdiff(reg_regions, lab_reg$region[lab_reg$variable == "Sea.Level"]), collapse = ", "),
    nrow(lab_noaa))
say("NOAA zones: %d U.S. regional zones (%s) plus the Mauna Loa CO2 series (World)",
    length(setdiff(unique(lab_noaa$region_raw), "World")), paste(setdiff(unique(lab_noaa$region_raw), "World"), collapse = "; "))
say("Climate components standardised on %s; macro transformations: %s", paste(range(eval(formals(standardize_monthly)$ref_years)), collapse = "-"),
    with(read.csv(file.path(PATHS$data, "macro_transformations.csv"), stringsAsFactors = FALSE),
         paste(sprintf("%s %s", series, transformation_used), collapse = ", ")))
r <- est$comp_national
say("Calendar %s to %s: %d months; COVID months %s to %s removed (%d) -> %d; estimation sample ends %s (%d months)",
    ym(min(r$dates_full)), ym(max(r$dates_full)), length(r$dates_full),
    ym(min(CFG$covid_dates)), ym(max(CFG$covid_dates)), length(CFG$covid_dates),
    length(r$dates), ym(r$end_test_date), sum(as.Date(r$dates) <= as.Date(r$end_test_date)))
say("VAR(%d) with constant, climate first: %d observations (supplied band objects: n_obs = %s), %d fitted residual rows, of which %d gap-dummied",
    r$lags, refit$comp_national$n_y, paste(unique(vapply(bands, `[[`, numeric(1), "n_obs")), collapse = "/"),
    refit$comp_national$n_resid, refit$comp_national$n_gap)
say("Penalty selected by %d-fold contiguous block CV; lambda-hat: %s", EST$cv_folds,
    paste(sprintf("%s %.6g", names(est), vapply(est, function(e) as.numeric(e$lambda), numeric(1))), collapse = ", "))
say("Bands: %.0f%% level", 100 * CFG$ci_level)

## ---- impulse responses ----------------------------------------------------------------------
block("[Sec. 3.2.1, Figs. 1-2; Sec. 1; Sec. 4] Trough of the cumulated INDPRO response (fixed-weights bands)")
b0 <- bands$benchmark$bands$INDPRO
for (id in names(bands)) {
  b <- bands[[id]]$bands$INDPRO; h <- which.min(b$irf)
  rf <- refit[[id]]$irf[, "INDPRO"]
  say("%-14s %.3f%% at h = %d (month %d), 68%% band [%.3f, %.3f], ratio to benchmark %.2f; band below zero at h = %s; re-fit from step 1: %.3f%% at h = %d (max |diff| %.1e)",
      id, b$irf[h], h - 1L, h, b$lower[h], b$upper[h], min(b$irf) / min(b0$irf), hrange(which(b$upper < 0) - 1L),
      min(rf), which.min(rf) - 1L, max(abs(refit[[id]]$irf - bands[[id]]$irf_point)))
}
tr <- vapply(bands, function(x) min(x$bands$INDPRO$irf), numeric(1))
say("comp_regional and comp_noaa troughs: %.2f%% and %.2f%% (%.3f / %.3f); %.2f and %.2f times the benchmark trough (%.2f%%)",
    -tr["comp_regional"], -tr["comp_noaa"], -tr["comp_regional"], -tr["comp_noaa"],
    tr["comp_regional"] / tr["benchmark"], tr["comp_noaa"] / tr["benchmark"], -tr["benchmark"])
say("Largest assembled trough relative to the benchmark: %.2f (%s)", max(tr[-1] / tr["benchmark"]), names(which.max(tr[-1] / tr["benchmark"])))
rr <- bands$rank_regional$bands$INDPRO$irf[1:12]; bb <- b0$irf[1:12]
say("rank_regional vs benchmark, months 1-12 (h = 0..11): mean %.3f%% vs %.3f%%, trough %.3f%% vs %.3f%%; rank_regional below the benchmark at h = %s",
    mean(rr), mean(bb), min(rr), min(bb), hrange(which(rr < bb) - 1L))

block("[Sec. 3.2.1, Figs. 1-2] Unemployment and housing starts: horizons h at which the 68% band excludes zero")
for (id in names(bands)) {
  u <- bands[[id]]$bands$UNRATE; hs <- bands[[id]]$bands$HOUST
  say("%-14s UNRATE above zero h = %s, below zero h = %s, point at h = 36 %+.3f | HOUST below zero h = %s, above zero h = %s; HOUST point > 0 at %d of 37 horizons",
      id, hrange(which(u$lower > 0) - 1L), hrange(which(u$upper < 0) - 1L), u$irf[37],
      hrange(which(hs$upper < 0) - 1L), hrange(which(hs$lower > 0) - 1L), sum(hs$irf > 0))
}

## ---- component weights ----------------------------------------------------------------------
block("[Sec. 1; Sec. 3.2.2 and Fig. 3, lower left] ACI component weights (point estimates)")
w <- est$comp_national$weights
say("National: %s (equal weight %s%%)", shares(w), pct(1 / length(w), 0))
w <- est$comp_regional$weights
say("Regional, by variable family: %s", shares(tapply(w, lab_reg$variable, sum)))
say("Regional, WP90 by region (share of the whole index): %s", shares(setNames(w[lab_reg$variable == "WP90"], lab_reg$region[lab_reg$variable == "WP90"])))
say("Regional, T90 by region: %s", shares(setNames(w[lab_reg$variable == "T90"], lab_reg$region[lab_reg$variable == "T90"])))
say("Regional, Rx5Day by region: %s", shares(setNames(w[lab_reg$variable == "Rx5Day"], lab_reg$region[lab_reg$variable == "Rx5Day"])))
say("Regional, total weight by region: %s", shares(tapply(w, lab_reg$region, sum)))
wind_share <- sapply(split(seq_along(w), lab_reg$region), function(i) sum(w[i][lab_reg$variable[i] == "WP90"]) / sum(w[i]))
say("Regional, wind share within each region's own weight: %s", shares(wind_share, 0))

block("[Sec. 3.2.2 and Fig. 4, lower left; App. C] NOAA component weights by region x variable cell (point estimates)")
w <- est$comp_noaa$weights
cell <- tapply(w, list(lab_noaa$region_agg, lab_noaa$variable_agg), sum); cell[is.na(cell)] <- 0
say("Cells (percent):"); tab <- capture.output(print(round(100 * cell, 1))); cat(tab, sep = "\n"); out <- c(out, tab)
say("Variable families: %s", shares(colSums(cell)))
say("Ohio Valley thermal block (Temp. + Heat.DD): %s%%; Northern Rockies + Northeast hydroclimate (Precip. + Palmer.Z): %s%%",
    pct(cell["Ohio.V.", "Temp."] + cell["Ohio.V.", "Heat.DD"], 2),
    pct(sum(cell[c("N.Rockies", "N.East"), c("Precip.", "Palmer.Z")]), 2))

## ---- rank weights ---------------------------------------------------------------------------
block("[Sec. 3.2.2 and Figs. 3-4, lower right; App. B, NOAA panel] Rank weights (point estimates)")
say("Tails = ranks below the 20th or above the 80th percentile; maxima of the raw profile and of the displayed (smoothed) profile, left and right halves")
for (id in rank_ids) {
  w <- abs(est[[id]]$weights) / sum(abs(est[[id]]$weights)); K <- length(w); q <- rank_percentile(K)
  sm <- smooth_rank_profile(w, id); n <- length(sm); qs <- plotted_percentile(n)
  L <- seq_len(floor(K / 2)); R <- setdiff(seq_len(K), L); Ls <- seq_len(floor(n / 2)); Rs <- setdiff(seq_len(n), Ls)
  jl <- L[which.max(w[L])]; jr <- R[which.max(w[R])]; il <- Ls[which.max(sm[Ls])]; ir <- Rs[which.max(sm[Rs])]
  say("%s: K = %d, flat weight %s%%; raw max left %s%% (rank %d, pct %.0f), right %s%% (rank %d, pct %.0f); displayed max left %s%% (pct %.0f), right %s%% (pct %.0f); largest single rank %s%%",
      id, K, pct(1 / K, 2), pct(w[jl], 3), jl, 100 * q[jl], pct(w[jr], 3), jr, 100 * q[jr],
      pct(sm[il], 3), 100 * qs[il], pct(sm[ir], 3), 100 * qs[ir], pct(max(w), 3))
  say("%s: tail mass %s%%, 30th-55th percentile mass %s%% (%s%% with percentile j/K), 25th-55th %s%%",
      id, pct(sum(w[q < 0.2 | q > 0.8]), 1), pct(rank_mass(w, 0.30, 0.55), 2),
      pct(sum(w[seq_len(K) / K >= 0.30 & seq_len(K) / K <= 0.55]), 2), pct(rank_mass(w, 0.25, 0.55), 2))
  cum <- vapply(seq(0.1, 0.9, 0.1), function(p) sum(w[q <= p]), numeric(1))
  say("%s: cumulative mass up to the 10th, 20th, ..., 90th percentile: %s", id, paste(pct(cum, 1), collapse = ", "))
}
say("rank_national raw profile (ranks 1-5): %s", paste(pct(est$rank_national$weights, 2), collapse = ", "))

## ---- Appendix A.1 ---------------------------------------------------------------------------
block("[App. A.1] Fixed-weights residual bootstrap (supplied band objects)")
for (id in names(bands)) {
  b <- bands[[id]]
  say("%-14s n_obs %d, fitted residuals %d (floor(n^(1/3)) = %d), block length %d, B_outer %d, B_inner %d, B_se %d, seed %d, all bands finite %s",
      id, b$n_obs, refit[[id]]$n_resid, floor(refit[[id]]$n_resid^(1 / 3)), b$block_length, b$B_outer, b$B_inner, b$B_se, b$seed,
      all(vapply(b$bands, function(z) all(is.finite(c(z$lower, z$upper))), logical(1))))
}

## ---- Appendix A.2 ---------------------------------------------------------------------------
block("[App. A, App. A.2 and notes to Figs. A.1-A.2] Weight-reestimating pairs bootstrap: draws")
for (id in names(SPECS)) {
  t1 <- task1[[id]]
  usable <- sum(t1$boot$converged %in% TRUE & !vapply(t1$boot$var_B, is.null, logical(1)))
  say("%-14s draws %d, weight solutions %d, usable VAR fits %d, pairs %d, block length %d = max(p + 1, floor(%d^(1/3)) = %d)",
      id, t1$n_draws, t1$n_converged, usable, t1$sample_n, t1$block_length, t1$sample_n, floor(t1$sample_n^(1 / 3)))
}

block("[App. A.2, Fig. A.1] Bootstrap weights: point, median, 68% and 90% intervals")
for (id in names(SPECS)) {
  ws <- bootstrap_weight_statistics(id, est[[id]]$weights, task1[[id]], lab_reg, lab_noaa)
  if (startsWith(id, "comp_")) {
    say("%s: %s", id, paste(sprintf("%s %s (med %s; 68%% %s-%s; 90%% %s-%s)", ws$group, pct(ws$point), pct(ws$q50),
                                    pct(ws$q16), pct(ws$q84), pct(ws$q05), pct(ws$q95)), collapse = "; "))
  } else {
    i <- which.max(ws$point); n <- nrow(ws)
    say("%s: displayed profile peaks at position %d of %d (rank j/K = %.0f%%, plotted pct %.0f) with %s%% (68%% %s-%s)",
        id, i, n, 100 * i / n, 100 * plotted_percentile(n)[i], pct(ws$point[i], 2), pct(ws$q16[i], 2), pct(ws$q84[i], 2))
  }
}
W <- usable_boot_weights(task1$comp_regional)
say("comp_regional: WP90.MID 68%% %s-%s, WP90.SPL 68%% %s-%s; WP90.MID > WP90.SPL in %s%% of draws",
    pct(quantile(W[, lab_reg$name == "WP90.MID"], 0.16)), pct(quantile(W[, lab_reg$name == "WP90.MID"], 0.84)),
    pct(quantile(W[, lab_reg$name == "WP90.SPL"], 0.16)), pct(quantile(W[, lab_reg$name == "WP90.SPL"], 0.84)),
    pct(mean(W[, lab_reg$name == "WP90.MID"] > W[, lab_reg$name == "WP90.SPL"]), 0))

block("[App. A.2, Fig. A.2] Point index inside its pointwise bootstrap bands (share of plotted dates)")
cov <- t(vapply(names(SPECS), function(id)
  index_band_coverage(bootstrap_index_paths(id, SPECS[[id]], PATHS$data, est[[id]]$weights, task1[[id]])), numeric(3)))
for (id in names(SPECS))
  say("%-14s %d dates: inside 68%% band %s%%, inside 90%% band %s%%", id, cov[id, "n"], pct(cov[id, "cover68"]), pct(cov[id, "cover90"]))
say("Average over component space: 68%% %s%%; over rank space: 68%% %s%%",
    pct(mean(cov[comp_ids, "cover68"])), pct(mean(cov[rank_ids, "cover68"])))

block("[App. A.2, Figs. A.3-A.4] INDPRO under the weight-reestimating pairs bootstrap (68% basic bands)")
for (id in names(SPECS)) {
  bo <- bands[[id]]
  cumulate <- vapply(bo$bands, function(b) isTRUE(b$cumulated), logical(1))
  pb <- pairs_bootstrap_irf(task1[[id]], var_order = bo$var_order, p = CFG$lags, horizon = CFG$irf_horizon,
                            cumulate = cumulate, macro_order = CFG$macro_order)
  ci <- pairs_bootstrap_bands(pb$point, pb$draws, conf = CFG$ci_level, reverse_percentile = TRUE)
  b <- bo$bands$INDPRO; h <- which.min(b$irf)
  say("%-14s %d draws; trough %.3f%% at h = %d, pairs band [%.3f, %.3f] (fixed-weights [%.3f, %.3f]); pairs band below zero at h = %s",
      id, pb$n_draws, b$irf[h], h - 1L, ci$lower["INDPRO", h], ci$upper["INDPRO", h], b$lower[h], b$upper[h],
      hrange(which(ci$upper["INDPRO", ] < 0) - 1L))
}

## ---- Appendix A.3 ---------------------------------------------------------------------------
block("[App. A.3, Fig. A.5] Expanding-window shares (COVID months at 2019-12 values, centred 12-month mean as plotted)")
for (id in names(SPECS)) {
  Wm <- task2[[id]]$W_mat
  if (startsWith(id, "comp_")) {
    groups <- appendix_groups(id, est[[id]]$weights, lab_reg, lab_noaa)
    lev <- if (id == "comp_national") names(est[[id]]$weights) else unique(groups)
  } else {
    Wm <- t(apply(Wm, 1L, smooth_rank_profile, id = id))   # quintiles of the displayed profile
    groups <- rank_quintile_groups(ncol(Wm)); lev <- paste0("Q", 1:5)
  }
  sh <- expanding_share_matrix(Wm, groups, lev); dt <- rownames(sh)
  say("%s (%s to %s): %s", id, ym(dt[1]), ym(dt[length(dt)]),
      paste(sprintf("%s %s->%s%% (mean %s, min %s %s, max %s %s)", colnames(sh), pct(sh[1, ]), pct(sh[nrow(sh), ]), pct(colMeans(sh)),
                    pct(apply(sh, 2, min)), ym(dt[apply(sh, 2, which.min)]),
                    pct(apply(sh, 2, max)), ym(dt[apply(sh, 2, which.max)])), collapse = "; "))
  if (id == "comp_national") {
    i <- substr(dt, 1, 4) >= "2010" & substr(dt, 1, 4) <= "2015"
    say("comp_national, 2010-2015: Rx5Day min %s%% (%s), WP90 min %s%%", pct(min(sh[i, "Rx5Day"])),
        ym(dt[i][which.min(sh[i, "Rx5Day"])]), pct(min(sh[i, "WP90"])))
  }
  if (id == "rank_regional") {
    yr <- substr(dt, 1, 4)
    for (g in c("Q4", "Q5"))
      say("rank_regional %s, annual means %s-%s: %s", g, yr[1], yr[length(yr)], paste(pct(tapply(sh[, g], yr, mean), 0), collapse = " "))
  }
}

## ---- Appendix B -----------------------------------------------------------------------------
block("[App. B, Fig. B.1] Rank-implied usage U of named components")
u <- rank_usage_noaa(est$rank_noaa$weights, est$comp_noaa$weights, est$rank_noaa$comps_perm, lab_noaa)
say("NOAA family U: %s; Spearman with the component-space family weights %.2f",
    paste(sprintf("%s %+.1f%%", u$family$family, 100 * u$family$usage), collapse = ", "), u$spearman_family)
u <- rank_usage_components(est$rank_national$weights, est$rank_national$comps_perm)
say("ACI national U by component: %s", paste(sprintf("%s %+.1f%%", names(est$comp_national$weights), 100 * u), collapse = ", "))
u <- rank_usage_components(est$rank_regional$weights, est$rank_regional$comps_perm)
say("ACI regional U, family means: %s; components range %+.1f%% to %+.1f%%",
    paste(sprintf("%s %+.1f%%", names(tapply(u, lab_reg$variable, mean)), 100 * tapply(u, lab_reg$variable, mean)), collapse = ", "),
    100 * min(u), 100 * max(u))

## ---- Appendix C -----------------------------------------------------------------------------
block("[App. C, Tables C.1-C.2] NOAA aggregation tables vs data/labels_noaa.csv")
say("Regions kept separately: %s; pooled into Rest: %s", paste(setdiff(unique(lab_noaa$region_agg), "Rest"), collapse = ", "),
    paste(unique(lab_noaa$region_raw[lab_noaa$region_agg == "Rest"]), collapse = "; "))
say("Variable families: %s; pooled into Others: %s", paste(setdiff(unique(lab_noaa$variable_agg), "Others"), collapse = ", "),
    paste(unique(lab_noaa$variable_raw[lab_noaa$variable_agg == "Others"]), collapse = "; "))
say("Number of raw NOAA series: %d; region x variable cells with at least one series: %d",
    nrow(lab_noaa), nrow(unique(lab_noaa[, c("region_agg", "variable_agg")])))

writeLines(out, file.path(PATHS$summary, "summary_numbers.txt"))
cat("\nwritten results/summary_numbers.txt\n")
