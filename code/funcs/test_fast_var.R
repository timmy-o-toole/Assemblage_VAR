## The direct least-squares VAR bootstrap refit must match vars::VAR numerically.
source("config.R")
source(file.path(ROOT, "code", "funcs", "irf_bands.R"))
suppressWarnings(suppressPackageStartupMessages(library(vars)))

cum <- c(climate = TRUE, INDPRO = TRUE, PCEPI = TRUE, UNRATE = TRUE,
         SHADOWRATE = FALSE, EXUSUKx = TRUE, HOUST = TRUE)
r <- readRDS(file.path(PATHS$estimates, "comp_national.rds"))
Y <- cbind(climate = r$index, as.matrix(r$macro[, CFG$macro_order]))
rownames(Y) <- r$dates

Yfit <- Y[2:(nrow(Y) - 1), ]
mv <- vars::VAR(Yfit, p = 12, type = "const")
mf <- var_fit_fast(Yfit, 12)
cv <- var_coefs(mv)
cf <- var_coefs(mf)
differences <- c(
  coefficients = max(abs(cv$A - cf$A)),
  constant = max(abs(cv$const - cf$const)),
  residuals = max(abs(stats::residuals(mv) - mf$resid)),
  roots = max(abs(sort(var_roots(mv)) - sort(var_roots(mf)))),
  irf = max(abs(irf_matrix(mv, "climate", 36L, cum) - irf_matrix(mf, "climate", 36L, cum)))
)

model <- fit_var(Y, 12L, r$dates, holdout = 1L)
resid <- scale(residuals(model), center = TRUE, scale = FALSE)
ok <- gap_ok_resid_rows(nrow(resid), 12L, gap_bad_positions(r$dates, 12L))
set.seed(123)
bootstrap_worst <- 0
for (b in 1:5) {
  nr <- mbb_resample(resid, 8, ok)
  mf <- refit_var_mbb(model, nr)
  cf <- var_coefs(model)
  K <- model$K; p <- model$p; y <- model$y; Tn <- nrow(y)
  yb <- matrix(NA_real_, Tn, K, dimnames = list(NULL, colnames(y)))
  yb[1:p, ] <- y[1:p, ]
  for (t in (p + 1):Tn) {
    x_lag <- as.vector(t(yb[(t - 1):(t - p), , drop = FALSE]))
    yb[t, ] <- cf$const + as.numeric(cf$A %*% x_lag) + nr[t - p, ]
  }
  mv <- vars::VAR(yb, p = p, type = "const")
  bootstrap_worst <- max(bootstrap_worst,
    max(abs(irf_matrix(mv, "climate", 36L, cum) - irf_matrix(mf, "climate", 36L, cum))),
    max(abs(stats::residuals(mv) - mf$resid)))
}
stopifnot(max(differences) < 1e-9, bootstrap_worst < 1e-9)
cat(sprintf("FAST VAR TEST PASSED: maximum difference %.2e\n", max(differences, bootstrap_worst)))
