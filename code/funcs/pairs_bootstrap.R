## =============================================================================
## code/funcs/pairs_bootstrap.R -- weight-reestimating pairs bootstrap (Task 1) and
## expanding-window re-estimation (Task 2) of the aggregation weights.
##
## Numerical core of 03_bootstrap_weights.R; it performs the computation behind the
## supplied results in results/bootstrap/. Wherever the package already has a
## numerically identical helper (lag_stack, eye_blocks, gap_valid_pairs, loglik_gain,
## build_design, penalty in code/funcs/estimation.R) it is reused; the sampler, the
## VAR fit and the origin bookkeeping are specific to this step.
##
## Conventions:
##  * Origins. The driver works on the lag stack trimmed by the one-month holdout
##    (T_full = est$n_use = 574 rows). Out-of-sample origin o = 1..216 uses rows
##    1..t_end with t_end = T_full - (217 - o); origin 216 is the fullest sample,
##    t_end = 573. The criterion is evaluated on the aligned (Y, X) pairs
##    r = 1..t_end-1 (pair r joins lag-stack rows r and r+1), restricted to the
##    pairs whose 13-month window does not cross the removed covid months. At
##    origin 216 this leaves 560 pairs.
##  * Constraints. w >= 0 and 1'w = 1 everywhere. Rank space adds the equal-mean
##    restriction sum(w * colMeans(comps)) = benchmark mean, where the benchmark
##    mean is taken over the whole trimmed sample for the ACI specifications
##    (est$bench_full) and over the lag-1 benchmark column of the rows available
##    at the origin for NOAA (est$bench_lag1[1:t_end]).
##  * Solver. NlcOptim::solnl with the settings EST in config.R; the
##    point estimate starts at the equal-weight vector, every bootstrap draw
##    starts at the point estimate of the same origin; a draw whose solve errors
##    is recorded as not converged.
##  * RNG. RNGkind("L'Ecuyer-CMRG"), set.seed(SEED) at the top of the script; the
##    block indices of origin o are drawn after set.seed(SEED + o) with the
##    global RNG state saved before and restored afterwards.
##  * Labels. The VAR matrices are labelled as follows: the
##    climate (assembled) variable is unnamed (""), regressor rows are
##    "(const)" then, per lag, the six macro variables followed by the assembled
##    lag; residual rows are labelled with the month of lag-stack row r (one month
##    before the dependent observation of pair r).
## =============================================================================

## ---- block length and gap-aware moving-block sampler ---------------------------------------

## Block length: never shorter than the VAR memory p.
choose_block_length <- function(T, p) {
  as.integer(max(p + 1L, floor(T^(1/3))))
}

## GAP-AWARE moving-block index on the PAIR space. `ok_pairs` is the ascending set of (Y, X)
## pairs whose lag window does not cross a removed-date hole. Blocks are drawn INSIDE contiguous
## runs of ok_pairs (runs weighted by length), so no resampled block straddles the hole and no
## excluded pair can re-enter through the bootstrap. Returns length(ok_pairs) pair indices.
make_mbb_pair_index_valid <- function(ok_pairs, ell) {
  runs <- unname(split(ok_pairs, cumsum(c(1L, diff(ok_pairs) != 1L))))
  len  <- vapply(runs, length, 1L)
  if (any(len < ell))
    warning(sprintf(paste0("make_mbb_pair_index_valid: %d run(s) shorter than block length %d ",
                           "cannot be sampled -> %d of %d valid pairs (%.1f%%) excluded from the ",
                           "bootstrap. Shorten the block or widen the sample."),
                    sum(len < ell), ell, sum(len[len < ell]), length(ok_pairs),
                    100 * sum(len[len < ell]) / length(ok_pairs)))
  runs <- runs[len >= ell]; len <- len[len >= ell]
  if (!length(runs))
    stop("make_mbb_pair_index_valid: no run of valid pairs >= block length ", ell)
  n <- length(ok_pairs); out <- integer(0)
  while (length(out) < n) {
    k <- sample.int(length(runs), 1L, prob = len)
    r <- runs[[k]]
    s <- if (length(r) == ell) 1L else sample.int(length(r) - ell + 1L, 1L)
    out <- c(out, r[s:(s + ell - 1L)])
  }
  out[seq_len(n)]
}

make_mbb_pair_indices_valid <- function(ok_pairs, B, ell) {
  lapply(seq_len(B), function(i) make_mbb_pair_index_valid(ok_pairs, ell))
}

## ---- origin bookkeeping -------------------------------------------------------------------------

## Rows and gap-valid pairs available at out-of-sample origin `origin` (1..n_origins).
origin_pairs <- function(est, origin, n_origins = 216L) {
  t_end   <- est$n_use - (n_origins + 1L - origin)         # last lag-stack row used
  n_pairs <- t_end - 1L                                    # aligned (Y, X) rows
  ok      <- intersect(seq_len(n_pairs), est$valid_pairs)  # minus pairs crossing the hole
  list(origin = origin, t_end = t_end, n_pairs = n_pairs, ok_pairs = ok)
}

## The estimation object restricted to the first t_end lag-stack rows.
truncate_est <- function(est, t_end) {
  est$macro_lags <- est$macro_lags[seq_len(t_end), , drop = FALSE]
  est$comps_lags <- est$comps_lags[seq_len(t_end), , drop = FALSE]
  est$t_end <- t_end
  est
}

## Target of the rank-space equal-mean constraint at an origin (see header).
constraint_bench_mean <- function(est, spec, t_end) {
  if (spec$level == 3L) mean(est$bench_lag1[seq_len(t_end)]) else mean(est$bench_full)
}

## ---- constrained solve on a set of pairs -------------------------------------------------------

## Minimise -S_T(w) + lambda * pen(w) over the pairs `rows` of the truncated design, starting at
## x0. Same objective as solve_weights() in code/estimation.R (loglik_gain + penalty) and the same
## solnl settings (tolX.cool2, tolFun 1e-2, tolCon 1e-8, 4000 iterations,
## 60000 evaluations); the start value and the benchmark mean are inputs here because the
## bootstrap needs them per draw / per origin. Returns the weight vector, or NULL when
## stop_on_error = FALSE and the solver errors (a non-converged bootstrap draw).
solve_weights_pairs <- function(est_o, spec, lambda, EST, rows, x0, bench_mean, stop_on_error = TRUE) {
  rank <- isTRUE(spec$rank); is_noaa <- spec$level == 3L
  K <- est_o$K; shrink_to <- rep(1 / K, K)
  objective <- function(w) {
    des <- build_design(est_o, w)
    loglik_gain(des$Z, des$take_col, ridge = EST$ridge_xtx, rows = rows) +
      penalty(w, lambda, rank, est_o, shrink_to, EST$penalty_scale)
  }
  Aeq <- matrix(1, 1, K); Beq <- 1
  eq_mean <- function(w) list(ceq = sum(w * est_o$comps_mean) - bench_mean, c = NULL)
  no_con  <- function(w) list(ceq = NULL, c = NULL)
  args <- list(X = x0, objfun = objective, lb = rep(0, K), Aeq = Aeq, Beq = Beq,
               maxnFun = EST$max_fun, maxIter = EST$max_iter,
               tolX = EST$tol_x_final(rank), tolFun = EST$tol_fun_final)
  if (rank) {                     # rank space: equal-mean constraint
    args <- c(args, list(confun = eq_mean, tolCon = EST$tol_con))
  } else if (is_noaa) {           # NOAA component space: constrained call with an empty constraint set
    args <- c(args, list(confun = no_con, tolCon = EST$tol_con))
  }                               # ACI component space: linear constraints only
  fit <- if (stop_on_error) do.call(NlcOptim::solnl, args)
         else tryCatch(do.call(NlcOptim::solnl, args), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  as.numeric(fit$par)
}

## ---- VAR on the aligned pairs ------------------------------------------------------------------

fit_var_sample <- function(Y_full, X_full, rows, ridge = 0, include_const = TRUE) {
  X_tr <- X_full[rows, , drop = FALSE]
  Y_tr <- Y_full[rows, , drop = FALSE]
  if (include_const) X_tr <- cbind(`(const)` = 1, X_tr)   # unpenalised constant regressor
  XtX <- crossprod(X_tr)
  XtY <- crossprod(X_tr, Y_tr)
  ridge_mat <- ridge * diag(ncol(X_tr))
  if (include_const) ridge_mat[1, 1] <- 0
  coeff <- tryCatch(
    solve(XtX + ridge_mat, XtY, tol = 1e-12),
    error = function(e) tryCatch(qr.solve(XtX + ridge_mat, XtY, tol = 1e-12), error = function(e2) NULL))
  if (is.null(coeff)) return(NULL)
  rownames(coeff) <- colnames(X_tr)
  colnames(coeff) <- colnames(Y_tr)
  residuals <- Y_tr - X_tr %*% coeff
  colnames(residuals) <- colnames(Y_tr)
  list(coeff = coeff, residuals = residuals, mu = colMeans(Y_tr))
}

## VAR(p) with constant for weights w on the pairs `rows` of the truncated sample. Regressor
## order per lag j: the six macro lags, then the assembled lag (85 rows incl. the constant);
## equations: the six macro variables, then the assembled variable. `dates` = d$dates.
fit_var_pairs <- function(est_o, w, rows, dates, ridge = 0) {
  p <- est_o$p; n_macro <- ncol(est_o$macro_lags); M <- n_macro / p
  assembled <- est_o$comps_lags %*% (t(est_o$eye) * rep(w, p))       # p lags of w'X, unnamed
  Z <- cbind(est_o$macro_lags, assembled)
  rownames(Z) <- dates[seq_len(nrow(Z)) + p - 1L]                    # row labels (see header)
  take_col  <- c(seq_len(M), n_macro + 1L)                            # lag-1 block + assembled lag 1
  take_col2 <- unlist(lapply(seq_len(p), function(j) c(((j - 1) * M + 1):(j * M), n_macro + j)))
  Y_full <- Z[2:nrow(Z), take_col, drop = FALSE]
  X_full <- Z[1:(nrow(Z) - 1), take_col2, drop = FALSE]
  fit_var_sample(Y_full, X_full, rows, ridge)
}

## ---- Task 1: point estimate and pairs bootstrap at one origin ---------------------------------

## Returns weights, VAR coefficients/residuals of the point sample and, for B > 0, the bootstrap
## draws (weights B x K, converged, corner, var_B and var_resid lists of length B; NULL entries
## for draws that did not converge).
bootstrap_origin <- function(est, d, spec, lambda, EST, origin, B, n_origins = 216L, seed = 123L,
                             verbose = TRUE) {
  op    <- origin_pairs(est, origin, n_origins)
  est_o <- truncate_est(est, op$t_end)
  bench_mean <- constraint_bench_mean(est, spec, op$t_end)
  K <- est$K; shrink_to <- rep(1 / K, K)
  series <- if (spec$level == 3L) "MOJO" else "ACI"
  if (verbose)
    cat(sprintf("  origin %d: rows 1..%d, %d pairs, %d gap-valid (%d cross the hole)\n",
                origin, op$t_end, op$n_pairs, length(op$ok_pairs), op$n_pairs - length(op$ok_pairs)))

  ## block indices for this origin: reproducible substream, global RNG state left untouched
  block_length <- NA_integer_; boot_pairs <- NULL
  if (B > 0) {
    rng_keep <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed + origin)
    block_length <- choose_block_length(length(op$ok_pairs), est$p)
    boot_pairs   <- make_mbb_pair_indices_valid(op$ok_pairs, B, block_length)
    assign(".Random.seed", rng_keep, envir = .GlobalEnv)
    if (verbose) cat(sprintf("  %d bootstrap draws, block length %d; first block starts at pairs %s\n",
                             B, block_length, paste(head(boot_pairs[[1]][c(1, block_length + 1)], 2), collapse = ", ")))
  }

  ## point sample: started at the equal-weight vector
  t0 <- Sys.time()
  w_point <- solve_weights_pairs(est_o, spec, lambda, EST, rows = op$ok_pairs, x0 = shrink_to,
                                 bench_mean = bench_mean, stop_on_error = TRUE)
  var_point <- fit_var_pairs(est_o, w_point, op$ok_pairs, d$dates, EST$ridge_xtx)
  if (is.null(var_point)) stop("Point VAR estimation failed for origin ", origin)
  if (verbose) cat(sprintf("  point weights solved in %.1f s; sum(w) = %.8f; largest: %s\n",
                           as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(w_point),
                           paste(sprintf("%.4f", head(sort(w_point, decreasing = TRUE), 5)), collapse = " ")))

  out <- list(origin = origin, t_end = op$t_end, sample_n = length(op$ok_pairs), block_length = block_length,
              weights = w_point, var_B = var_point$coeff, var_resid = var_point$residuals, boot = NULL)

  ## bootstrap draws: started at the point estimate, same penalty and constraints
  W <- matrix(NA_real_, B, K, dimnames = list(NULL, rep(series, K)))
  converged <- rep(FALSE, B); corner <- rep(NA, B)
  B_list <- vector("list", B); R_list <- vector("list", B)
  t1 <- Sys.time()
  for (b in seq_len(B)) {
    wb <- solve_weights_pairs(est_o, spec, lambda, EST, rows = boot_pairs[[b]], x0 = w_point,
                              bench_mean = bench_mean, stop_on_error = FALSE)
    if (!is.null(wb)) {
      W[b, ] <- wb; converged[b] <- TRUE
      corner[b] <- any(abs(wb) < 1e-08)
      vb <- fit_var_pairs(est_o, wb, boot_pairs[[b]], d$dates, EST$ridge_xtx)
      if (!is.null(vb)) { B_list[[b]] <- vb$coeff; R_list[[b]] <- vb$residuals }
    }
    if (verbose && (b %% 10 == 0 || b == B)) {
      el <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
      cat(sprintf("  draw %d/%d  converged %d  elapsed %.0f s  eta %.0f min\n",
                  b, B, sum(converged[seq_len(b)]), el, el / b * (B - b) / 60))
    }
  }
  if (B > 0) out$boot <- list(weights = W, converged = converged, corner = corner, var_B = B_list, var_resid = R_list)
  out
}

## ---- Task 2: expanding-window point weights -------------------------------------------------

## Point weights re-estimated at each origin in `origins` (subset of 1..n_origins), each solve
## started at the equal-weight vector. Returns the length(origins) x K matrix W_mat.
expanding_weights <- function(est, d, spec, lambda, EST, origins, n_origins = 216L, verbose = TRUE) {
  K <- est$K; shrink_to <- rep(1 / K, K)
  series <- if (spec$level == 3L) "MOJO" else "ACI"
  W <- matrix(NA_real_, length(origins), K, dimnames = list(as.character(origins), rep(series, K)))
  t0 <- Sys.time()
  for (i in seq_along(origins)) {
    o <- origins[i]
    op <- origin_pairs(est, o, n_origins)
    est_o <- truncate_est(est, op$t_end)
    W[i, ] <- solve_weights_pairs(est_o, spec, lambda, EST, rows = op$ok_pairs, x0 = shrink_to,
                                  bench_mean = constraint_bench_mean(est, spec, op$t_end), stop_on_error = TRUE)
    if (verbose) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("  origin %3d (rows 1..%d, %d pairs): sum(w) = %.6f, %d weights > 1e-6  [%.0f s, eta %.0f min]\n",
                  o, op$t_end, length(op$ok_pairs), sum(W[i, ]), sum(W[i, ] > 1e-6), el,
                  el / i * (length(origins) - i) / 60))
    }
  }
  W
}
