## =============================================================================
## code/funcs/cv.R -- non-overlapping block cross-validation of the penalty lambda.
##
## The gap-valid (Y, X) pairs are cut into `folds` contiguous blocks. For each
## lambda on the grid and each fold m, the weights are estimated on the other
## folds (training pairs must be adjacent, so pairs straddling a fold boundary or
## the covid hole are dropped) and the held-out criterion is the log-likelihood
## gain of a VAR estimated on the held-out block alone. The selected lambda
## minimises the fold average of Log0 - Log1 (i.e. maximises the average gain).
## Deterministic: no random numbers are used, so the result does not depend on
## the number of workers.
## =============================================================================

## Held-out criterion on pairs `pairs_out` of the full design Z (VAR refitted on the block).
holdout_gain <- function(Z, take_col, pairs_out, ridge = 0) {
  Zo   <- Z[pairs_out, , drop = FALSE]
  keep <- which(diff(pairs_out) == 1L)                 # adjacent original rows only
  if (!length(keep)) stop("no adjacent hold-out pairs in this fold")
  Y <- as.matrix(Zo[2:nrow(Zo), take_col, drop = FALSE])[keep, , drop = FALSE]
  X <- as.matrix(Zo[1:(nrow(Zo) - 1), , drop = FALSE])[keep, , drop = FALSE]
  XtX <- crossprod(X); XtY <- crossprod(X, Y)
  b <- solve(XtX + ridge * diag(ncol(X)), XtY, tol = 1e-12)
  residual <- as.matrix(Y)
  for (j in 1:ncol(Y)) residual[, j] <- residual[, j] - mean(residual[, j])
  Sig0 <- crossprod(residual) / nrow(residual); Inv0 <- solve(Sig0, tol = 1e-12)
  Log0 <- - nrow(residual) * log(det(Sig0))
  for (t in 1:nrow(residual)) Log0 <- Log0 - residual[t, ] %*% Inv0 %*% residual[t, ]
  residual <- as.matrix(Y - X %*% b)
  Sig1 <- crossprod(residual) / nrow(residual); Inv1 <- solve(Sig1, tol = 1e-12)
  Log1 <- - nrow(residual) * log(det(Sig1))
  for (t in 1:nrow(residual)) Log1 <- Log1 - residual[t, ] %*% Inv1 %*% residual[t, ]
  as.numeric(Log0 - Log1)
}

## One (lambda, fold) cell: estimate on the training pairs, score on the held-out block.
cv_cell <- function(est, d, spec, lambda, fold, folds_id, pool, EST) {
  d_in  <- pool[folds_id != fold]
  d_out <- pool[folds_id == fold]
  keep_in <- which(diff(d_in) == 1L)
  if (!length(keep_in)) stop("no adjacent training pairs in this fold")
  ## training: the criterion is evaluated on Z[d_in, ] restricted to adjacent pairs
  est_in <- est
  est_in$macro_lags <- est$macro_lags; est_in$comps_lags <- est$comps_lags
  rows_fun <- function(Z) Z[d_in, , drop = FALSE]
  fit <- solve_weights_rows(est, d, spec, lambda, EST, d_in, keep_in)
  des <- build_design(est, fit$weights)
  list(value = holdout_gain(des$Z, des$take_col, d_out, EST$ridge_xtx), weights = fit$weights)
}

## Variant of solve_weights() whose criterion is evaluated on the sub-design Z[d_in, ] with
## pairs `keep_in` (training folds). Uses the CV solver tolerances and the lag-1 benchmark mean.
solve_weights_rows <- function(est, d, spec, lambda, EST, d_in, keep_in) {
  rank <- isTRUE(spec$rank); is_noaa <- spec$level == 3L
  K <- est$K; shrink_to <- rep(1 / K, K)
  objective <- function(w) {
    des <- build_design(est, w)
    s <- loglik_gain(des$Z[d_in, , drop = FALSE], des$take_col, ridge = EST$ridge_xtx, rows = keep_in)
    s + penalty(w, lambda, rank, est, shrink_to, EST$penalty_scale)
  }
  Aeq <- matrix(1, 1, K); Beq <- 1
  bench_mean <- mean(est$bench_lag1)
  eq_mean <- function(w) list(ceq = sum(w * est$comps_mean) - bench_mean, c = NULL)
  no_con  <- function(w) list(ceq = NULL, c = NULL)
  common <- list(X = shrink_to, objfun = objective, lb = rep(0, K), Aeq = Aeq, Beq = Beq,
                 maxnFun = EST$max_fun, maxIter = EST$max_iter, tolX = EST$tol_x_cv(rank))
  args <- if (rank) c(common, list(confun = eq_mean, tolCon = EST$tol_con))
          else if (is_noaa) c(common, list(confun = no_con, tolCon = EST$tol_con))
          else common
  fit <- do.call(NlcOptim::solnl, args)
  list(weights = as.numeric(fit$par))
}

## Full grid x folds cross-validation. Returns the selected lambda, its grid index, and the
## fold-score matrix (grid x folds).
cross_validate <- function(est, d, spec, grid, EST, n_cores = 1L, verbose = TRUE) {
  ## As in the paper, the CV uses the full calendar, including the final month
  ## that is held out from the point estimate.
  est <- prepare_estimation(d, est$p, holdout = 0L)
  pool <- sort(est$valid_pairs)
  nf <- EST$cv_folds
  fd <- rep(1:(nf - 1), each = floor(length(pool) / nf))
  fd <- c(fd, rep(nf, length(pool) - length(fd)))
  cells <- expand.grid(g = seq_along(grid), f = seq_len(nf))
  one <- function(i) {
    value <- tryCatch({
      r <- cv_cell(est, d, spec, grid[cells$g[i]], cells$f[i], fd, pool, EST)
      if (length(r$value) != 1L || !is.finite(r$value)) stop("non-finite CV score")
      as.numeric(r$value)
    }, error = function(e) NA_real_)
    c(g = cells$g[i], f = cells$f[i], value = value)
  }
  if (verbose) cat(sprintf("  CV: %d lambdas x %d folds = %d fits (%d workers)\n", length(grid), nf, nrow(cells), n_cores))
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("est", "d", "spec", "grid", "cells", "fd", "pool", "EST",
                                  "cv_cell", "solve_weights_rows", "holdout_gain", "build_design",
                                  "loglik_gain", "penalty"), envir = environment())
    res <- parallel::parLapply(cl, seq_len(nrow(cells)), one)
  } else {
    res <- lapply(seq_len(nrow(cells)), one)
  }
  res <- do.call(rbind, res)
  scores <- matrix(NA_real_, length(grid), nf)
  scores[cbind(res[, "g"], res[, "f"])] <- res[, "value"]
  complete <- rowSums(is.finite(scores)) == nf
  fold_mean <- rep(NA_real_, length(grid))
  fold_mean[complete] <- rowMeans(scores[complete, , drop = FALSE])
  n_failed_cells <- sum(!is.finite(scores))
  n_failed_lambdas <- sum(!complete)
  if (n_failed_cells > 0L)
    warning(sprintf("CV discarded %d failed cell(s) and %d incomplete lambda row(s)",
                    n_failed_cells, n_failed_lambdas), call. = FALSE)
  if (!any(complete)) stop("CV failed: no lambda completed all folds", call. = FALSE)
  idx <- which.min(replace(fold_mean, !complete, Inf))
  list(lambda = grid[idx], grid_index = idx, fold_scores = scores, fold_mean = fold_mean, grid = grid,
       fold_sizes = as.integer(table(fd)), n_pool = length(pool))
}
