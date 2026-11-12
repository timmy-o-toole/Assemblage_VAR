## =============================================================================
## 02_irf_bands.R -- Step 2: impulse responses and fixed-weights
## residual bootstrap bands for the six specifications and for the equal-weight
## benchmark (Section 3, Appendix A.1).
##
## Reads results/point_estimates/<spec>.rds. With RUN_IRF_BOOTSTRAP = FALSE
## (default), the supplied band objects are checked and retained. With
## RUN_IRF_BOOTSTRAP = TRUE, the full nested bootstrap writes only to
## results/irf_bands/_recomputed/ and never replaces supplied results
## (1000 x 100 + 200 VAR refits per VAR; about 1 h each on one core).
##
## Usage:  Rscript code/02_irf_bands.R [spec_id ... | benchmark]
##         RUN_IRF_BOOTSTRAP=TRUE Rscript code/02_irf_bands.R
##         B_OUTER=20 B_INNER=5 B_SE=5 RUN_IRF_BOOTSTRAP=TRUE ... (smoke test)
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "irf_bands.R"))
suppressWarnings(suppressPackageStartupMessages(library(vars)))

args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args)) args else c(names(SPECS), "benchmark")
B_OUTER <- as.integer(Sys.getenv("B_OUTER", "1000")); B_INNER <- as.integer(Sys.getenv("B_INNER", "100")); B_SE <- as.integer(Sys.getenv("B_SE", "200"))

cumulate <- c(climate = TRUE, INDPRO = TRUE, PCEPI = TRUE, UNRATE = TRUE, SHADOWRATE = FALSE, EXUSUKx = TRUE, HOUST = TRUE)

## VAR data for one specification: [assembled index | macro variables] on the estimation sample
var_data <- function(r, benchmark = FALSE) {
  x <- if (benchmark) r$benchmark else r$index
  Y <- cbind(climate = x, as.matrix(r$macro[, CFG$macro_order]))
  rownames(Y) <- r$dates
  Y
}

for (id in ids) {
  is_bench <- id == "benchmark"
  ref_id <- if (is_bench) "comp_national" else id      # the benchmark uses the national sample (identical dates for all specs)
  r <- readRDS(file.path(PATHS$estimates, paste0(ref_id, ".rds")))
  Y <- var_data(r, benchmark = is_bench)
  holdout <- length(r$dates) - match(r$end_test_date, r$dates)
  cat(sprintf("\n== %s ==  VAR(%d) on %d observations (rows 2..%d), climate first\n", id, r$lags, nrow(Y) - holdout - 1, nrow(Y) - holdout))
  model <- fit_var(Y, r$lags, r$dates, holdout = holdout, use_dummies = !is_bench)
  bad_pos <- gap_bad_positions(r$dates, r$lags)
  ok_resid <- gap_ok_resid_rows(nrow(residuals(model)), r$lags, bad_pos, first_row = 2L)
  irf_point <- irf_matrix(model, "climate", CFG$irf_horizon, cumulate)
  out_file <- if (RUN_IRF_BOOTSTRAP) {
    dir.create(PATHS$irf_recomputed, recursive = TRUE, showWarnings = FALSE)
    file.path(PATHS$irf_recomputed, paste0(id, ".rds"))
  } else {
    irf_result_file(id)
  }

  if (RUN_IRF_BOOTSTRAP) {
    set.seed(SEED)
    t0 <- Sys.time()
    bands <- percentile_t_bands(model, "climate", CFG$irf_horizon, cumulate, ok_resid,
                                B_outer = B_OUTER, B_inner = B_INNER, B_se = B_SE, conf = CFG$ci_level,
                                inner_se_floor = !is_bench)
    cat(sprintf("  done in %.1f min; rejected draws: %d; outer draws missing: %d\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), attr(bands, "rejected_draws"), attr(bands, "n_missing_outer")))
    cat("  redraws after a non-stationary refit (outer, inner, s.e.):", paste(attr(bands, "redraws"), collapse = ", "), "\n")
    saveRDS(list(id = id, bands = bands, irf_point = irf_point, var_order = colnames(Y),
                 n_obs = nrow(model$y), block_length = attr(bands, "block_length"),
                 B_outer = B_OUTER, B_inner = B_INNER, B_se = B_SE, seed = SEED,
                 rejected_draws = attr(bands, "rejected_draws"), redraws = attr(bands, "redraws"), session = list(time = Sys.time(), R = R.version.string, vars = as.character(packageVersion("vars")))),
            out_file)
    cat("  saved", sub(ROOT, ".", out_file, fixed = TRUE), "\n")
  } else {
    require_result(out_file, paste0("IRF band result for ", id))
    supplied <- readRDS(out_file)
    if (is.null(supplied$bands) || !identical(dim(supplied$irf_point), dim(irf_point)))
      stop("invalid IRF band result: ", out_file, call. = FALSE)
    cat("  using ", result_source(out_file), " residual-bootstrap bands\n", sep = "")
  }
}
