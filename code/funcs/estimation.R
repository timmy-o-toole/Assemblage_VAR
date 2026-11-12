## =============================================================================
## code/funcs/estimation.R -- Assemblage VAR: build the lagged design, evaluate the
## learning criterion, and solve for the aggregation weights.
##
## Notation follows the paper: the assembled variable x*_t(w) = w' X_t replaces
## the climate benchmark in a VAR(p) with constant; w minimises
##    -S_T(w) + lambda * pen(w)   s.t.  w >= 0, 1'w = 1,
## where S_T is the in-sample Gaussian quasi-log-likelihood gain of the VAR over
## a mean-only model, evaluated on (Y, X) pairs whose lag window does not cross
## the removed covid months.
## =============================================================================

## ---- lag stacking ---------------------------------------------------------------------------

## Rows (p+1)..n of [X_{t-1} | X_{t-2} | ... | X_{t-p}]  (n-p rows, K*p columns).
lag_stack <- function(x, p) {
  x <- as.matrix(x); n <- nrow(x)
  out <- do.call(cbind, lapply(seq_len(p), function(j) x[(p + 1 - j):(n - j), , drop = FALSE]))
  rownames(out) <- rownames(x)[(p + 1):n]
  out
}

## Selector matrix E (p x K*p): row j has ones on the columns of lag block j.
## t(E) * w (columnwise recycling) gives the K*p x p matrix whose product with the
## lag stack returns the p lags of the assembled series w'X.
eye_blocks <- function(K, p) {
  E <- matrix(0, p, K * p)
  for (j in seq_len(p)) E[j, ((j - 1) * K + 1):(j * K)] <- 1
  E
}

## ---- calendar gaps ----------------------------------------------------------------------------

## Pair r of the lag stack joins stack rows r (regressors) and r+1 (dependent). It is valid
## only if the p+1 calendar months it spans are consecutive. Returns the valid pair indices.
gap_valid_pairs <- function(dates, p) {
  per  <- 12L * as.integer(format(as.Date(dates), "%Y")) + as.integer(format(as.Date(dates), "%m"))
  d    <- diff(per)
  step <- as.integer(names(which.max(table(d))))
  brk  <- as.integer(d != step)
  cum  <- c(0L, cumsum(brk))
  n_rows <- length(dates) - p
  pair <- seq_len(n_rows - 1L)
  pair[cum[pair + p] == cum[pair]]
}

## ---- learning criterion -------------------------------------------------------------------

## Gaussian quasi-log-likelihood of a mean-only model minus that of the one-step VAR
## projection, both evaluated at their OLS estimates on the pairs `rows` of Z.
## Z: design rows t = 1..n (lag stack + assembled lags + constant); take_col: columns of Z whose
## next-row values are the dependent variables. Returns Log0 - Log1 (to be minimised).
loglik_gain <- function(Z, take_col, ridge = 0, rows = NULL, tol_fit = 1e-12, tol_sig = 1e-12) {
  Yraw <- as.matrix(Z[2:nrow(Z), take_col, drop = FALSE])
  X    <- as.matrix(Z[1:(nrow(Z) - 1), , drop = FALSE])
  if (!is.null(rows)) { Yraw <- Yraw[rows, , drop = FALSE]; X <- X[rows, , drop = FALSE] }
  Y <- Yraw
  XtX   <- crossprod(X)
  XtY   <- crossprod(X, Y)
  betas <- solve(XtX + ridge * diag(ncol(X)), XtY, tol = tol_fit)
  residual <- Y
  for (j in 1:ncol(Y)) residual[, j] <- residual[, j] - mean(residual[, j])
  Sig0 <- crossprod(residual) / nrow(residual)
  InvSig0 <- tryCatch(solve(Sig0, tol = tol_sig), error = function(e) qr.solve(Sig0, tol = tol_sig))
  Log0 <- - nrow(residual) * log(det(Sig0))
  for (t in 1:nrow(residual)) Log0 <- Log0 - residual[t, ] %*% InvSig0 %*% residual[t, ]
  residual <- Y - X %*% betas
  Sig1 <- crossprod(residual) / nrow(residual)
  InvSig1 <- tryCatch(solve(Sig1, tol = tol_sig), error = function(e) qr.solve(Sig1, tol = tol_sig))
  Log1 <- - nrow(residual) * log(det(Sig1))
  for (t in 1:nrow(residual)) Log1 <- Log1 - residual[t, ] %*% InvSig1 %*% residual[t, ]
  as.numeric(Log0 - Log1)
}

## ---- estimation sample -----------------------------------------------------------------------

## Build everything the optimiser needs from a loaded specification (see code/funcs/load_data.R).
## The last `holdout` months are excluded from the estimation window (holdout = 1
## in the paper, i.e. the sample ends 2025:M4).
prepare_estimation <- function(d, p, holdout = NULL) {
  if (is.null(holdout)) holdout <- length(d$dates) - match(d$end_test_date, d$dates)
  macro_lags <- lag_stack(d$macro, p)                 # (n-p) x (6p)
  bench_lags <- lag_stack(d$bench, p)                 # (n-p) x (7p), untrimmed
  comps_lags <- lag_stack(d$comps, p)                 # (n-p) x (Kp)
  n_use <- nrow(macro_lags) - holdout
  valid <- gap_valid_pairs(d$dates, p)
  list(
    macro_lags = macro_lags[1:n_use, , drop = FALSE],
    comps_lags = comps_lags[1:n_use, , drop = FALSE],
    bench_lag1 = bench_lags[, d$climate_name],       # lag-1 column of the climate benchmark (untrimmed)
    bench_full = d$bench[[d$climate_name]],           # benchmark over the trimmed sample
    comps_mean = colMeans(d$comps),
    comps_sd   = apply(d$comps, 2, stats::sd),
    K = ncol(d$comps), p = p, n_use = n_use,
    take_col   = 1:ncol(d$macro),
    eye        = eye_blocks(ncol(d$comps), p),
    valid_pairs_upto = function(n_rows) valid[valid <= (n_rows - 1L)],
    valid_pairs = valid,
    holdout = holdout
  )
}

## Design matrix for weights w: [macro lags | p lags of w'X | 1]; dependent columns = take_col + assembled.
build_design <- function(est, w) {
  assembled_lags <- est$comps_lags %*% (t(est$eye) * rep(w, est$p))   # weights recycled over the p lag blocks
  Z <- cbind(est$macro_lags, assembled_lags, 1)
  list(Z = Z, take_col = c(est$take_col, ncol(est$macro_lags) + 1L))
}

## Penalty: component space shrinks towards the equal-weight reference (scaled by the component
## standard deviations); rank space penalises differences between adjacent rank weights.
penalty <- function(w, lambda, rank, est, shrink_to, scale = 1) {
  if (rank) lambda * scale * sum(diff(est$comps_sd * w)^2)
  else      lambda * scale * sum(est$comps_sd * (w - shrink_to)^2)
}

## ---- solver ---------------------------------------------------------------------------------

## Solve for the weights of one specification given lambda.
##   spec$rank, spec$level; `mode` = "final" or "cv" selects the solver tolerances and, for the
##   rank-space equal-mean constraint, the benchmark mean:
##   final ACI: mean over the trimmed sample; CV ACI and NOAA (both modes): mean of the lag-1 column.
##   `rows` restricts the criterion to given pairs (used by the CV); default = gap-valid pairs.
solve_weights <- function(est, d, spec, lambda, EST, mode = c("final", "cv"), rows = NULL, verbose = FALSE) {
  mode <- match.arg(mode)
  rank <- isTRUE(spec$rank); is_noaa <- spec$level == 3L
  K <- est$K
  shrink_to <- rep(1 / K, K)
  if (is.null(rows)) rows <- est$valid_pairs_upto(nrow(est$macro_lags))
  n_eval <- 0L
  objective <- function(w) {
    n_eval <<- n_eval + 1L
    des <- build_design(est, w)
    s <- loglik_gain(des$Z, des$take_col, ridge = EST$ridge_xtx, rows = rows)
    if (verbose && n_eval %% 500 == 0) cat(sprintf("    %d evaluations, criterion %.3f\n", n_eval, s))
    s + penalty(w, lambda, rank, est, shrink_to, EST$penalty_scale)
  }
  ## equality constraints: 1'w = 1 (linear); rank space adds the equal-mean restriction
  Aeq <- matrix(1, 1, K); Beq <- 1
  bench_mean <- if (rank && !is_noaa && mode == "final") mean(est$bench_full) else mean(est$bench_lag1)
  eq_mean <- function(w) list(ceq = sum(w * est$comps_mean) - bench_mean, c = NULL)
  no_con  <- function(w) list(ceq = NULL, c = NULL)
  lb <- rep(0, K); x0 <- shrink_to
  tolX <- if (mode == "final") EST$tol_x_final(rank) else EST$tol_x_cv(rank)
  common <- list(X = x0, objfun = objective, lb = lb, Aeq = Aeq, Beq = Beq,
                 maxnFun = EST$max_fun, maxIter = EST$max_iter, tolX = tolX)
  if (mode == "final") common$tolFun <- EST$tol_fun_final
  if (rank) {                     # rank space: equal-mean constraint, constrained SQP
    args <- c(common, list(confun = eq_mean, tolCon = EST$tol_con))
  } else if (is_noaa) {           # NOAA component space: constrained call with an empty constraint set
    args <- c(common, list(confun = no_con, tolCon = EST$tol_con))
  } else {                        # ACI component space: linear constraints only
    args <- common
  }
  fit <- do.call(NlcOptim::solnl, args)
  w <- as.numeric(fit$par); names(w) <- colnames(d$comps)
  list(weights = w, objective = as.numeric(fit$fn), n_eval = n_eval, lambda = lambda, mode = mode)
}

## Assembled index over the trimmed sample and the criterion at the solution.
assembled_index <- function(d, w) {
  x <- as.numeric(as.matrix(d$comps) %*% w); names(x) <- rownames(d$comps); x
}
