## =============================================================================
## config.R -- all settings of the empirical application in one place.
## Sourced by run_all.R and by the numbered scripts in code/.
## =============================================================================

## ---- run switches (can be overridden from the command line / run_all.R) --------------------
## All four default to FALSE: the fast route reruns no expensive estimator and regenerates
## none of the archived pairs-bootstrap and expanding-window estimation outputs.
env_flag <- function(name, default = "FALSE") {
  v <- toupper(trimws(Sys.getenv(name, default)))
  if (v %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO", "")) return(FALSE)
  stop("environment variable ", name, " must be TRUE or FALSE, got '", Sys.getenv(name), "'", call. = FALSE)
}
RUN_CV                <- env_flag("RUN_CV")                # re-select the penalty by block CV
RUN_IRF_BOOTSTRAP     <- env_flag("RUN_IRF_BOOTSTRAP")     # recompute the fixed-weights residual bootstrap bands
RUN_BOOTSTRAP_WEIGHTS <- env_flag("RUN_BOOTSTRAP_WEIGHTS") # re-run the appendix pairs bootstrap / expanding windows
USE_RECOMPUTED        <- env_flag("USE_RECOMPUTED")        # read fresh IRF and appendix results only
N_CORES           <- as.integer(Sys.getenv("N_CORES", max(1L, parallel::detectCores() - 1L)))
SEED              <- 123L

## ---- sample and VAR -------------------------------------------------------------------------
CFG <- list(
  macro_order   = c("INDPRO", "PCEPI", "UNRATE", "SHADOWRATE", "EXUSUKx", "HOUST"), # Cholesky order after the climate variable
  start_date    = "1975-01-01",
  end_date      = "2025-05-01",          # last month in the data; the estimation window ends one month earlier
  end_test_date = "2025-04-01",          # last month used in estimation (one trailing month is held out)
  covid_dates   = as.character(seq.Date(as.Date("2020-01-01"), as.Date("2021-06-01"), by = "1 month")),
  lags          = 12L,
  irf_horizon   = 36L,
  ci_level      = 0.68
)

## ---- the six specifications ----------------------------------------------------------------
## level: 1 = ACI national (5 series), 2 = ACI regional (34), 3 = NOAA (82)
## rank : FALSE = component space, TRUE = rank space
## grid_index: position of the cross-validated penalty of the published estimates on the CV grid
## (block CV, 5 folds). The penalty itself is the exact grid value, filled in below; it is used when
## RUN_CV = FALSE. Values: 11.56838, 5.056576, 26.46603, 7.389056, 49.23239, 141785.3.
SPECS <- list(
  comp_national = list(id = "comp_national", level = 1L, rank = FALSE, panel = "ACI national", label = "National, component space",  grid_index = 8L),
  comp_regional = list(id = "comp_regional", level = 2L, rank = FALSE, panel = "ACI regional", label = "Regional, component space",  grid_index = 4L),
  comp_noaa     = list(id = "comp_noaa",     level = 3L, rank = FALSE, panel = "NOAA",         label = "NOAA, component space",      grid_index = 12L),
  rank_national = list(id = "rank_national", level = 1L, rank = TRUE,  panel = "ACI national", label = "National, rank space",       grid_index = 1L),
  rank_regional = list(id = "rank_regional", level = 2L, rank = TRUE,  panel = "ACI regional", label = "Regional, rank space",       grid_index = 6L),
  rank_noaa     = list(id = "rank_noaa",     level = 3L, rank = TRUE,  panel = "NOAA",         label = "NOAA, rank space",           grid_index = 27L)
)

## ---- estimator settings ---------------------------------------------------------------------
EST <- list(
  ridge_xtx     = 0,                       # ridge added to X'X in the VAR score (none)
  cv_folds      = 5L,                      # contiguous blocks
  penalty_scale = 1,                       # multiplier on the penalty
  ## CV grid for the penalty: component space and rank space use different ranges
  lambda_grid   = function(rank) if (rank) exp(seq(2, 13, length.out = 30)) else exp(seq(1, 7, length.out = 30)),
  ## SQP (NlcOptim::solnl) settings, as used for the published estimates
  max_iter      = 4000L,
  max_fun       = 60000L,
  tol_x_cv      = function(rank) if (rank) 1e-05 else 1e-04,   # tolX inside the CV
  tol_x_final   = function(rank) if (rank) 1e-08 else 1e-06,   # tolX for the final solve
  tol_fun_final = 1e-02,
  tol_con       = 1e-08
)
for (nm in names(SPECS)) SPECS[[nm]]$lambda <- EST$lambda_grid(SPECS[[nm]]$rank)[SPECS[[nm]]$grid_index]

## ---- paths -----------------------------------------------------------------------------------
find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = TRUE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "config.R")) && dir.exists(file.path(d, "code"))) return(d)
    d <- dirname(d)
  }
  stop("package root (folder containing config.R and code/) not found above ", start)
}
ROOT <- find_root()
PATHS <- list(
  data       = file.path(ROOT, "data"),
  estimates  = file.path(ROOT, "results", "point_estimates"),
  ## supplied, published residual-bootstrap results (read-only by convention)
  irf_bands  = file.path(ROOT, "results", "irf_bands"),
  ## fresh residual-bootstrap results; never replaces the supplied files
  irf_recomputed = file.path(ROOT, "results", "irf_bands", "_recomputed"),
  ## the archived pairs-bootstrap and expanding-window estimation outputs (read-only by convention)
  bootstrap  = file.path(ROOT, "results", "bootstrap"),
  ## a fresh run of 03_bootstrap_weights.R lands here and never touches the archive
  recomputed = file.path(ROOT, "results", "bootstrap", "_recomputed"),
  figures    = file.path(ROOT, "figures"),
  summary    = file.path(ROOT, "results")
)
## Fresh-result directories are created on demand by the estimators, not here,
## so a fresh package folder does not acquire empty result directories.
for (nm in setdiff(names(PATHS), c("recomputed", "irf_recomputed")))
  dir.create(PATHS[[nm]], recursive = TRUE, showWarnings = FALSE)

## ---- supplied vs recomputed outputs -----------------------------------------------------------
## Result selection is strict. If USE_RECOMPUTED is TRUE, callers receive a path in the relevant
## _recomputed directory even when the file is absent; downstream scripts then stop clearly.
boot_result_file <- function(id, task) {
  f <- paste0(id, "_task", task, ".rds")
  file.path(if (USE_RECOMPUTED) PATHS$recomputed else PATHS$bootstrap, f)
}

irf_result_file <- function(id) {
  file.path(if (USE_RECOMPUTED) PATHS$irf_recomputed else PATHS$irf_bands,
            paste0(id, ".rds"))
}

result_source <- function(path) {
  fresh_roots <- c(PATHS$recomputed, PATHS$irf_recomputed)
  if (any(startsWith(normalizePath(path, winslash = "/", mustWork = FALSE),
                     normalizePath(fresh_roots, winslash = "/", mustWork = FALSE))))
    "recomputed" else "supplied"
}
boot_result_source <- result_source

require_result <- function(path, description = "result") {
  if (!file.exists(path))
    stop(description, " not found: ", sub(ROOT, ".", path, fixed = TRUE),
         if (USE_RECOMPUTED) " (USE_RECOMPUTED=TRUE requires a complete fresh result set)" else "",
         call. = FALSE)
  path
}

point_estimate_lambda <- function(id, expected_k, estimate_dir = PATHS$estimates) {
  path <- require_result(file.path(estimate_dir, paste0(id, ".rds")),
                         paste0("point estimate for ", id))
  x <- readRDS(path)
  if (!identical(x$id, id) || length(x$weights) != expected_k ||
      length(x$lambda) != 1L || !is.finite(x$lambda) || x$lambda <= 0)
    stop("point-estimate metadata does not match specification ", id, call. = FALSE)
  as.numeric(x$lambda)
}
