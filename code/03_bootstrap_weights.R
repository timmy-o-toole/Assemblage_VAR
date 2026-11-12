## =============================================================================
## 03_bootstrap_weights.R -- Step 3: uncertainty and stability of the
## aggregation weights (appendix figures boot_weights_*, boot_index_*,
## expanding_weights_*).
##
## Task 1: at the final out-of-sample origin (216 = the full sample) estimate the
##         point weights with the lambda selected in step 1, then B weight-reestimating
##         pairs-bootstrap draws (moving blocks of gap-valid (Y, X) pairs); the
##         VAR is fitted on the same pairs for every draw.
## Task 2: re-estimate the point weights at each of the 216 expanding origins
##         (monthly from 2005-10 to 2025-03 on the covid-free calendar).
##
## The task-specific numerical routines are in code/funcs/pairs_bootstrap.R.
##
## Usage:  Rscript code/03_bootstrap_weights.R <spec_id> <task: 1|2> [B] [n_origins]
##           spec_id    one of names(SPECS) in config.R
##           B          bootstrap draws for task 1 (default 1000; the published run)
##           n_origins  task 2: run only the LAST n_origins of the 216 origins (default 216)
## Output: results/bootstrap/_recomputed/<spec_id>_task1.rds / _task2.rds
##         (format: results/README.md). A full run NEVER writes into results/bootstrap/
##         itself: those twelve files are the supplied pairs-bootstrap and expanding-window
##         results used in the paper and are treated as read-only.
##         Reduced runs (B < 1000 or n_origins < 216) are written to
##         results/bootstrap/_smoke/ so that they never overwrite a full result;
##         the environment variable BOOT_OUT_DIR overrides the output directory
##         (but may not point at the archive).
## Runs offline from the package root; one core is enough.
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "load_data.R"))
source(file.path(ROOT, "code", "funcs", "estimation.R"))
source(file.path(ROOT, "code", "funcs", "pairs_bootstrap.R"))

FINAL_ORIGIN <- 216L      # origin of the fullest sample
N_ORIGINS    <- 216L      # size of the out-of-sample origin grid
B_FULL       <- 1000L     # draws of the published run

args <- commandArgs(trailingOnly = TRUE)
usage <- "usage: Rscript code/03_bootstrap_weights.R <spec_id> <task: 1|2> [B] [n_origins]"
if (length(args) < 2) stop(usage, call. = FALSE)
id <- args[1]
if (!id %in% names(SPECS)) stop("unknown spec_id '", id, "'; use one of: ", paste(names(SPECS), collapse = ", "), call. = FALSE)
task <- as.integer(args[2])
if (!task %in% 1:2) stop(usage, call. = FALSE)
B         <- if (length(args) >= 3) as.integer(args[3]) else B_FULL
n_origins <- if (length(args) >= 4) as.integer(args[4]) else N_ORIGINS
if (task == 1L && (is.na(B) || B < 0)) stop("B must be a non-negative integer", call. = FALSE)
if (task == 2L && (is.na(n_origins) || n_origins < 1 || n_origins > N_ORIGINS)) stop("n_origins must be in 1..216", call. = FALSE)

full_run <- if (task == 1L) B == B_FULL else n_origins == N_ORIGINS
out_dir <- Sys.getenv("BOOT_OUT_DIR", "")
if (!nzchar(out_dir)) out_dir <- if (full_run) PATHS$recomputed else file.path(PATHS$bootstrap, "_smoke")
## The archived outputs in results/bootstrap/ are never a write target, whatever BOOT_OUT_DIR says.
norm_dir <- function(p) sub("/+$", "", normalizePath(p, winslash = "/", mustWork = FALSE))
if (identical(norm_dir(out_dir), norm_dir(PATHS$bootstrap)))
  stop("refusing to overwrite the archived pairs-bootstrap and expanding-window estimation ",
       "outputs in results/bootstrap/; a fresh run belongs in results/bootstrap/_recomputed/ ",
       "(unset BOOT_OUT_DIR, or point it somewhere else)", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("  output directory: %s\n", out_dir))

## RNG scheme: parallel-safe streams, seed 123 at the top; the block indices
## of origin o are drawn under set.seed(123 + o) inside bootstrap_origin().
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

spec <- SPECS[[id]]
cat(sprintf("\n== %s (%s) -- task %d%s ==\n", id, spec$label, task,
            if (task == 1L) sprintf(", B = %d", B) else sprintf(", last %d of %d origins", n_origins, N_ORIGINS)))
t_start <- Sys.time()
cfg <- CFG; cfg$climate_name <- "climate"
d   <- load_specification(spec, PATHS$data, cfg)
est <- prepare_estimation(d, CFG$lags)
lambda <- point_estimate_lambda(id, est$K)
cat(sprintf("  sample %s..%s, %d months without covid; lag stack %d rows after the %d-month holdout; K = %d; lambda = %.6g\n",
            min(d$dates), max(d$dates), length(d$dates), est$n_use, est$holdout, est$K, lambda))
session <- list(time = Sys.time(), R = R.version.string, platform = R.version$platform,
                NlcOptim = as.character(packageVersion("NlcOptim")), rng = RNGkind()[1], seed = SEED)

if (task == 1L) {
  res <- bootstrap_origin(est, d, spec, lambda, EST, origin = FINAL_ORIGIN, B = B, n_origins = N_ORIGINS, seed = SEED)
  K <- est$K
  boot <- if (B > 0) res$boot else
    list(weights = matrix(NA_real_, 0, K, dimnames = list(NULL, rep(if (spec$level == 3L) "MOJO" else "ACI", K))),
         converged = logical(0), corner = logical(0), var_B = list(), var_resid = list())
  out <- list(
    id = id, origin = FINAL_ORIGIN, lambda = lambda,
    lambda_source = file.path("results", "point_estimates", paste0(id, ".rds")),
    weights = res$weights,                                  # point weights at the final origin (560 pairs)
    var_B = res$var_B, var_resid = res$var_resid,           # VAR(12) with constant on the same pairs
    boot = list(weights = boot$weights, converged = boot$converged, var_B = boot$var_B, var_resid = boot$var_resid,
                corner = boot$corner),
    n_draws = B, n_converged = sum(boot$converged), sample_n = res$sample_n, block_length = res$block_length,
    note = sprintf("Weight-reestimating pairs bootstrap at the final origin (B = %d), 03_bootstrap_weights.R.", B),
    session = session
  )
  f <- file.path(out_dir, paste0(id, "_task1.rds"))
  saveRDS(out, f)
  cat(sprintf("  task 1 done: %d of %d draws converged, sample_n %d, block length %d; %.1f min; saved %s\n",
              out$n_converged, B, out$sample_n, out$block_length,
              as.numeric(difftime(Sys.time(), t_start, units = "mins")), f))
} else {
  origins <- seq.int(N_ORIGINS - n_origins + 1L, N_ORIGINS)
  W <- expanding_weights(est, d, spec, lambda, EST, origins = origins, n_origins = N_ORIGINS)
  out <- list(
    id = id, lambda = lambda,
    lambda_source = file.path("results", "point_estimates", paste0(id, ".rds")),
    origins = origins,
    W_mat = W,                                               # row o = weights with the sample ending at origin o
    weights_final = as.numeric(W[nrow(W), ]),
    note = sprintf("Expanding-window re-estimation of the point weights (%d origins), 03_bootstrap_weights.R.", length(origins)),
    session = session
  )
  f <- file.path(out_dir, paste0(id, "_task2.rds"))
  saveRDS(out, f)
  cat(sprintf("  task 2 done: %d origins (%d..%d); %.1f min; saved %s\n", length(origins), min(origins), max(origins),
              as.numeric(difftime(Sys.time(), t_start, units = "mins")), f))
}
