## =============================================================================
## code/funcs/irf_bands.R -- structural impulse responses and fixed-weights residual
## bootstrap bands (Appendix A.1 of the paper).
##
## Given the assembled climate index (weights held fixed), a VAR(p) with constant
## is fitted on the sample with the covid months removed; the p observations whose
## lag window crosses the hole are excluded through pulse dummies (algebraically
## identical to dropping them). Shocks are identified recursively (Cholesky) with
## the climate variable ordered first. Bands are pointwise equal-tailed
## percentile-t intervals from a moving-block residual bootstrap: B_outer
## artificial samples generated recursively from the fitted VAR, each studentised
## with B_inner inner draws; the standard error of the original response comes
## from B_se additional draws. Draws with a non-stationary refit are redrawn up to
## 50 times and otherwise recorded as missing.
##
## The original model is fitted with vars::VAR (with the gap dummies). The
## bootstrap refits use a direct least-squares fit on the same design matrix
## (var_fit_fast), which returns the same coefficients, residuals and impulse
## responses as vars::VAR/vars::irf up to floating-point rounding (checked in
## code/funcs/test_fast_var.R) and is about seven times faster. The sequence of random
## draws matches the supplied bands.
## =============================================================================

## ---- calendar helpers -----------------------------------------------------------------------

## Calendar positions (1..n) whose p-lag window crosses a hole in `dates`.
gap_bad_positions <- function(dates, p) {
  per  <- 12L * as.integer(format(as.Date(dates), "%Y")) + as.integer(format(as.Date(dates), "%m"))
  d    <- diff(per); step <- as.integer(names(which.max(table(d))))
  cum  <- c(0L, cumsum(as.integer(d != step)))
  obs  <- seq.int(p + 1L, length(dates))
  obs[cum[obs] != cum[obs - p]]
}

## Residual rows (of a VAR fitted on sample rows first_row:..., p lags) that are gap-clean.
gap_ok_resid_rows <- function(n_rows, p, bad_pos, first_row = 2L) {
  if (!length(bad_pos)) return(seq_len(n_rows))
  bad <- bad_pos - (first_row - 1L) - p
  setdiff(seq_len(n_rows), bad[bad >= 1L & bad <= n_rows])
}

## Admissible moving-block starts: the whole block must lie inside one gap-clean run.
mbb_starts <- function(n, ell, ok = NULL) {
  s <- seq_len(max(n - ell + 1L, 1L))
  if (is.null(ok) || length(ok) >= n) return(s)
  okv <- logical(n); okv[ok] <- TRUE
  keep <- vapply(s, function(z) all(okv[z:min(z + ell - 1L, n)]), logical(1))
  if (!any(keep)) stop("no block of length ", ell, " fits inside a gap-clean run")
  s[keep]
}

## ---- VAR fit -------------------------------------------------------------------------------

## Fit VAR(p) with constant on rows 2..(n - holdout) of Y (a T x K matrix, climate first),
## excluding the gap-contaminated observations through pulse dummies. Returns a vars::varest.
fit_var <- function(Y, p, dates, holdout = 1L, use_dummies = TRUE) {
  n <- nrow(Y); rows <- 2:(n - holdout)
  exogen <- NULL
  if (use_dummies) {
    pos <- gap_bad_positions(dates, p); pos <- pos[pos >= 1 & pos <= n]
    if (length(pos)) {
      D <- matrix(0, n, length(pos), dimnames = list(NULL, paste0("excl_", seq_along(pos))))
      for (j in seq_along(pos)) D[pos[j], j] <- 1
      E <- D[rows, , drop = FALSE]; E <- E[, colSums(E) > 0, drop = FALSE]
      if (ncol(E)) exogen <- E
    }
  }
  vars::VAR(Y[rows, , drop = FALSE], p = p, type = "const", exogen = exogen)
}

## Direct least-squares VAR(p) with constant (no exogenous terms): the same regression as
## vars::VAR(type = "const") -- lm.fit on [y_{t-1}, ..., y_{t-p}, 1] -- without the formula overhead.
var_fit_fast <- function(y, p) {
  y <- as.matrix(y); Tn <- nrow(y); K <- ncol(y)
  X <- do.call(cbind, lapply(seq_len(p), function(j) y[(p + 1 - j):(Tn - j), , drop = FALSE]))
  X <- cbind(X, 1)
  Yd <- y[(p + 1):Tn, , drop = FALSE]
  fit <- stats::lm.fit(X, Yd)
  coef <- fit$coefficients                          # (Kp+1) x K
  resid <- fit$residuals; dimnames(resid) <- list(NULL, colnames(y))
  structure(list(y = y, p = p, K = K, obs = Tn - p,
                 A = t(coef[seq_len(K * p), , drop = FALSE]),        # K x Kp, lag blocks in order
                 const = as.numeric(coef[K * p + 1, ]),
                 resid = resid), class = "var_fast")
}

## Accessors that work for both vars::varest and var_fast objects
var_y     <- function(m) if (inherits(m, "var_fast")) m$y else m$y
var_p     <- function(m) m$p
var_K     <- function(m) m$K
var_resid <- function(m) if (inherits(m, "var_fast")) m$resid else stats::residuals(m)
## K x Kp lag coefficient matrix and constant vector
var_coefs <- function(m) {
  if (inherits(m, "var_fast")) return(list(A = m$A, const = m$const))
  B <- vars::Bcoef(m)
  A <- B[, grep("\\.l[0-9]+$", colnames(B)), drop = FALSE]
  const_col <- which(colnames(B) == "const")
  list(A = A, const = if (length(const_col) == 1L) B[, const_col] else rep(0, m$K))
}
## Moduli of the companion-matrix eigenvalues (as vars::roots)
var_roots <- function(m) {
  if (!inherits(m, "var_fast")) return(Mod(vars::roots(m)))
  K <- m$K; p <- m$p; Fm <- matrix(0, K * p, K * p)
  Fm[1:K, ] <- m$A
  if (p > 1) Fm[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  Mod(eigen(Fm, only.values = TRUE)$values)
}

## ---- impulse responses --------------------------------------------------------------------

## Orthogonalised (Cholesky, in column order) responses of all variables to `impulse` over
## horizons 0..n_ahead, cumulated for the variables flagged in `cumulate` (named logical).
## For varest objects this is vars::irf; for var_fast objects the same computation is done directly
## (MA coefficients from the companion recursion; Sigma = U'U / (obs - Kp - 1) as in vars).
irf_matrix <- function(model, impulse, n_ahead, cumulate) {
  if (inherits(model, "var_fast")) {
    K <- model$K; p <- model$p
    A <- lapply(seq_len(p), function(j) model$A[, ((j - 1) * K + 1):(j * K), drop = FALSE])
    Sigma <- crossprod(model$resid) / (model$obs - K * p - 1)
    P <- t(chol(Sigma))
    j0 <- match(impulse, colnames(model$y))
    Phi <- vector("list", n_ahead + 1); Phi[[1]] <- diag(K)
    for (h in seq_len(n_ahead)) {
      Ph <- matrix(0, K, K)
      for (j in seq_len(min(h, p))) Ph <- Ph + A[[j]] %*% Phi[[h - j + 1]]
      Phi[[h + 1]] <- Ph
    }
    m <- t(vapply(Phi, function(Ph) (Ph %*% P)[, j0], numeric(K)))
    colnames(m) <- colnames(model$y)
  } else {
    o <- vars::irf(model, impulse = impulse, response = colnames(model$y), n.ahead = n_ahead,
                   boot = FALSE, cumulative = FALSE)
    m <- o$irf[[impulse]]
  }
  for (k in seq_len(ncol(m))) if (isTRUE(cumulate[colnames(m)[k]])) m[, k] <- cumsum(m[, k])
  m
}

## ---- bootstrap machinery (same random-draw sequence as the supplied bands) ------------------

mbb_resample <- function(resid_mat, block_len, ok_rows = NULL) {
  n <- nrow(resid_mat); n_blocks <- ceiling(n / block_len)
  starts <- sample(mbb_starts(n, block_len, ok_rows), size = n_blocks, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) s:(s + block_len - 1)))[1:n]
  resid_mat[idx, ]
}

## Generate an artificial sample recursively from the fitted VAR with resampled residuals,
## refit with the same lag order and constant, reject non-stationary refits (returns NULL).
refit_var_mbb <- function(model, new_resid, reject_unstable = TRUE) {
  K <- var_K(model); p <- var_p(model); y <- var_y(model); Tn <- nrow(y)
  cf <- var_coefs(model); A <- cf$A; c_hat <- cf$const
  y_boot <- matrix(NA_real_, Tn, K, dimnames = list(NULL, colnames(y)))
  y_boot[1:p, ] <- y[1:p, ]
  for (t in (p + 1):Tn) {
    x_lag <- as.vector(t(y_boot[(t - 1):(t - p), , drop = FALSE]))
    y_boot[t, ] <- c_hat + as.numeric(A %*% x_lag) + new_resid[t - p, ]
  }
  m <- var_fit_fast(y_boot, p)
  if (reject_unstable && any(var_roots(m) >= 1)) return(NULL)
  m
}

## Percentile-t bands. `ok_resid` restricts the outer block starts to gap-clean residual rows;
## the inner draws resample the outer residuals without restriction.
## `inner_se_floor`: TRUE = mask horizons with a numerically zero inner s.e. (assemblage VARs);
## FALSE = add 1e-10 to the s.e. (benchmark VAR). Both as for the supplied bands.
percentile_t_bands <- function(model, impulse, n_ahead, cumulate, ok_resid,
                               B_outer = 1000L, B_inner = 100L, B_se = 200L, conf = 0.68,
                               inner_se_floor = TRUE, max_attempts = 50L, verbose = TRUE) {
  alpha <- 1 - conf
  resid_orig <- scale(var_resid(model), center = TRUE, scale = FALSE)
  block_l <- max(floor(nrow(resid_orig)^(1/3)), 0)
  K <- var_K(model); resp <- colnames(var_y(model))
  irf_orig <- irf_matrix(model, impulse, n_ahead, cumulate)
  t_star <- replicate(K, matrix(NA_real_, B_outer, n_ahead + 1), simplify = FALSE); names(t_star) <- resp
  rejected <- 0L
  redraws <- c(outer = 0L, inner = 0L, se = 0L)   # extra draws after a non-stationary refit
  if (verbose) { cat(sprintf("  outer bootstrap: %d draws x %d inner, block length %d\n", B_outer, B_inner, block_l)); pb <- txtProgressBar(0, B_outer, style = 3) }
  for (b in seq_len(B_outer)) {
    m_out <- refit_var_mbb(model, mbb_resample(resid_orig, block_l, ok_resid))
    attempt <- 0L
    while (is.null(m_out) && attempt < max_attempts) { attempt <- attempt + 1L; m_out <- refit_var_mbb(model, mbb_resample(resid_orig, block_l, ok_resid)) }
    redraws[["outer"]] <- redraws[["outer"]] + attempt
    if (is.null(m_out)) { rejected <- rejected + 1L; next }
    irf_out <- irf_matrix(m_out, impulse, n_ahead, cumulate)
    resid_inner <- var_resid(m_out)
    irf_inner <- array(NA_real_, c(B_inner, n_ahead + 1, K))
    for (j in seq_len(B_inner)) {
      m_in <- refit_var_mbb(m_out, mbb_resample(resid_inner, block_l))
      attempt <- 0L
      while (is.null(m_in) && attempt < max_attempts) { attempt <- attempt + 1L; m_in <- refit_var_mbb(m_out, mbb_resample(resid_inner, block_l)) }
      redraws[["inner"]] <- redraws[["inner"]] + attempt
      if (is.null(m_in)) { rejected <- rejected + 1L; next }
      irf_inner[j, , ] <- irf_matrix(m_in, impulse, n_ahead, cumulate)
    }
    for (k in seq_len(K)) {
      se_inner <- apply(irf_inner[, , k], 2, stats::sd)
      if (inner_se_floor) {
        valid <- is.finite(se_inner) & se_inner > sqrt(.Machine$double.eps)
        t_star[[k]][b, valid]  <- (irf_out[valid, k] - irf_orig[valid, k]) / se_inner[valid]
        t_star[[k]][b, !valid] <- NA_real_
      } else {
        t_star[[k]][b, ] <- (irf_out[, k] - irf_orig[, k]) / (se_inner + 1e-10)
      }
    }
    if (verbose) setTxtProgressBar(pb, b)
  }
  if (verbose) { close(pb); cat("  standard error of the original response from", B_se, "refits\n") }
  se_mat <- do.call(rbind, lapply(seq_len(B_se), function(j) {
    m <- refit_var_mbb(model, mbb_resample(resid_orig, block_l, ok_resid))
    attempt <- 0L
    while (is.null(m) && attempt < max_attempts) { attempt <- attempt + 1L; m <- refit_var_mbb(model, mbb_resample(resid_orig, block_l, ok_resid)) }
    redraws[["se"]] <<- redraws[["se"]] + attempt
    if (is.null(m)) { rejected <<- rejected + 1L; return(rep(NA_real_, (n_ahead + 1) * K)) }
    as.vector(irf_matrix(m, impulse, n_ahead, cumulate))
  }))
  irf_se <- array(apply(se_mat, 2, stats::sd, na.rm = TRUE), dim = c(n_ahead + 1, K))
  out <- list()
  for (k in seq_len(K)) {
    q_lo <- apply(t_star[[k]], 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    q_hi <- apply(t_star[[k]], 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
    out[[resp[k]]] <- list(irf = irf_orig[, k], lower = irf_orig[, k] - q_lo * irf_se[, k],
                           upper = irf_orig[, k] - q_hi * irf_se[, k], cumulated = isTRUE(cumulate[resp[k]]))
  }
  attr(out, "n_missing_outer") <- sum(is.na(t_star[[1]][, 2]))
  attr(out, "rejected_draws")  <- rejected
  attr(out, "redraws")         <- redraws
  attr(out, "block_length")    <- block_l
  out
}
