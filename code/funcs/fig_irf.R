## =============================================================================
## code/funcs/fig_irf.R -- impulse-response panel figures (main text, Section 4, and
## Appendix) and the weight-reestimating pairs-bootstrap bands of the Appendix.
##
## One figure shows one specification: seven stacked panels, the climate index
## first and then the six macro variables in Cholesky order, each with the
## response to a one-standard-deviation climate shock over 36 months.
##   red  : the assembled index with its 68% band
##   blue : the equal-weight headline ACI benchmark with its 68% band
## The main-text figures (04a_figures_irf.R) draw the fixed-weights
## residual-bootstrap bands stored in results/irf_bands/<spec>.rds. The appendix
## figures (05a_figures_irfboot.R) replace the bands of the assembled
## index by pairs-bootstrap bands that also account for the estimation of the
## weights: the point response and the response of every converged bootstrap
## draw are recomputed here from the VAR coefficient matrices and residuals
## saved in results/bootstrap/<spec>_task1.rds.
##
## Requires: ggplot2, ggh4x (one y scale per panel), ggthemes (through
## code/funcs/theme.R). No fonts are downloaded; see irf_font() below.
## =============================================================================

## ---- labels and axis limits -----------------------------------------------------------------

## Strip (facet title) text on the left of each panel. The figures use the FRED-MD
## variable codes and show the climate index as "ACI"; VAR_LABELS in code/funcs/theme.R holds the long
## names used in captions and can be passed to draw_irf_panel(strip_labels = VAR_LABELS).
IRF_STRIP_LABELS <- c(climate = "ACI", INDPRO = "INDPRO", PCEPI = "PCEPI", UNRATE = "UNRATE",
                      SHADOWRATE = "SHADOWRATE", EXUSUKx = "EXUSUKx", HOUST = "HOUST")

## Short names printed inside the panels of the level-1 (national) figures, which stand at the
## left of a row of three figures in the paper (the regional and NOAA figures carry no labels).
IRF_PANEL_LABELS <- c(climate = "ACI", INDPRO = "Indust. Prod.", PCEPI = "Inflation", UNRATE = "Unemp.",
                      SHADOWRATE = "Shadow Rate", EXUSUKx = "Ex. Rate", HOUST = "Housing")

## Text scaling: draw_irf_panel() keeps the base size of the figure theme (45 pt) for margins and
## paddings and scales the explicit text sizes (axis text 50, strip text 55, ...) by this factor.
IRF_TEXT_SCALE <- 96 / 600

## Font family: "sans" resolves to Arial on Windows and
## macOS and to the default sans-serif face elsewhere, without any download.
irf_font <- function() "sans"

## y-axis limits and breaks per response variable. They differ between component space and
## rank space; the appendix (pairs-bootstrap) figures use slightly wider limits for three panels.
irf_y_scales <- function(rank, appendix = FALSE) {
  sc <- function(limits, breaks) ggplot2::scale_y_continuous(limits = limits, breaks = breaks)
  if (!rank) {
    list(
      climate    = sc(c(-0.05, 0.49), c(0, 0.1, 0.2, 0.3, 0.4)),
      INDPRO     = sc(c(-0.62, 0.14), c(-0.6, -0.45, -0.3, -0.15, 0)),
      PCEPI      = sc(c(-0.38, 0.35), c(-0.3, -0.15, 0, 0.15, 0.3)),
      UNRATE     = sc(c(-0.08, 0.15), c(-0.05, 0, 0.05, 0.1, 0.15)),
      SHADOWRATE = sc(if (appendix) c(-0.32, 0.32) else c(-0.3, 0.3), c(-0.3, -0.15, 0, 0.15, 0.3)),
      EXUSUKx    = sc(c(-0.6, 1.4), c(-0.4, 0, 0.4, 0.8, 1.2)),
      HOUST      = if (appendix) sc(c(-4.5, 1.8), c(-4.5, -3, -1.5, 0, 1.5)) else sc(c(-3.2, 1.5), c(-3, -2, -1, 0, 1))
    )
  } else {
    list(
      climate    = sc(c(-0.05, 0.49), c(0, 0.1, 0.2, 0.3, 0.4)),
      INDPRO     = sc(c(-0.48, if (appendix) 0.19 else 0.18), c(-0.4, -0.3, -0.2, -0.1, 0, 0.1)),
      PCEPI      = sc(c(-0.3, 0.15), c(-0.3, -0.2, -0.1, 0, 0.1)),
      UNRATE     = sc(c(-0.07, 0.08), c(-0.06, -0.03, 0, 0.03, 0.06)),
      SHADOWRATE = sc(c(-0.3, 0.18), c(-0.3, -0.2, -0.1, 0, 0.1)),
      EXUSUKx    = sc(c(-1, 1.1), c(-1, -0.5, 0, 0.5, 1)),
      HOUST      = sc(c(-1.5, 2.15), c(-1.5, -0.75, 0, 0.75, 1.5))
    )
  }
}

## In-panel labels (data coordinates) for a specification; NULL when the figure carries none.
irf_panel_labels <- function(level, rank, appendix = FALSE) {
  if (level != 1L) return(NULL)
  y <- if (!rank) {
    if (appendix) c(0.45, 0.11, 0.30, 0.13, 0.27, 1.18, 1.32) else c(0.45, 0.11, 0.30, 0.13, 0.25, 1.20, 1.35)
  } else {
    if (appendix) c(0.45, 0.15, 0.12, 0.07, 0.14, 0.90, 1.85) else c(0.45, 0.15, 0.12, 0.07, 0.14, 0.99, 1.85)
  }
  data.frame(response = names(IRF_PANEL_LABELS), label = unname(IRF_PANEL_LABELS),
             x = c(3.5, 10, 6.7, 6, 10, 7, 7.5), y = y, stringsAsFactors = FALSE)
}

## ---- plot data ------------------------------------------------------------------------------

## Long data frame (model, response, h, irf, lower, upper) from one band object
## (results/irf_bands/<id>.rds: `bands` = named list over the variables with irf/lower/upper).
irf_band_frame <- function(band_obj, model) {
  if (is.null(band_obj$bands)) stop("band object '", band_obj$id, "' has no bands (run 02_irf_bands.R with RUN_IRF_BOOTSTRAP = TRUE)")
  do.call(rbind, lapply(band_obj$var_order, function(v) {
    b <- band_obj$bands[[v]]
    data.frame(model = model, response = v, h = seq_along(b$irf) - 1L,
               irf = as.numeric(b$irf), lower = as.numeric(b$lower), upper = as.numeric(b$upper),
               stringsAsFactors = FALSE)
  }))
}

## Plot data for one figure: benchmark rows first, then the specification. `model` and
## `response` are factors (benchmark first; responses in Cholesky order).
irf_plot_data <- function(spec_bands, bench_bands, spec_name = "Assembled index", bench_name = "Headline ACI") {
  if (!identical(spec_bands$var_order, bench_bands$var_order)) stop("variable order differs between specification and benchmark")
  df <- rbind(irf_band_frame(bench_bands, bench_name), irf_band_frame(spec_bands, spec_name))
  df$model    <- factor(df$model, levels = c(bench_name, spec_name))
  df$response <- factor(df$response, levels = spec_bands$var_order)
  df
}

## ---- the figure -----------------------------------------------------------------------------

## df: from irf_plot_data(). y_scales: from irf_y_scales(). labels: from irf_panel_labels() or
## NULL. The layout reproduces the published 3 x 9.5 inch panel (aspect ratio 0.7 per panel).
draw_irf_panel <- function(df, y_scales, labels = NULL, strip_labels = IRF_STRIP_LABELS,
                           impulse_label = "ACI", base_family = irf_font(), text_scale = IRF_TEXT_SCALE) {
  models    <- levels(df$model)                      # benchmark, then the assembled index
  var_order <- levels(df$response)
  palette   <- stats::setNames(c("#004E89", "firebrick3"), models)
  df$impulse <- factor(impulse_label)
  bench <- df[df$model == models[1], ]
  spec  <- df[df$model == models[2], ]
  s <- text_scale

  p <- ggplot2::ggplot(df, ggplot2::aes(x = h, y = irf, group = model)) +
    ## benchmark: light ribbon, thin band edges, medium line
    ggplot2::geom_ribbon(data = bench, ggplot2::aes(ymin = lower, ymax = upper, fill = model), alpha = 0.06, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(data = bench, ggplot2::aes(y = lower, colour = model), linewidth = 0.08, alpha = 0.6, na.rm = TRUE) +
    ggplot2::geom_line(data = bench, ggplot2::aes(y = upper, colour = model), linewidth = 0.08, alpha = 0.6, na.rm = TRUE) +
    ggplot2::geom_line(data = bench, ggplot2::aes(colour = model), linewidth = 0.5, lineend = "round", na.rm = TRUE) +
    ## assembled index: stronger ribbon and edges, thick line
    ggplot2::geom_ribbon(data = spec, ggplot2::aes(ymin = lower, ymax = upper, fill = model), alpha = 0.10, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(data = spec, ggplot2::aes(y = lower, colour = model), linewidth = 0.12, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_line(data = spec, ggplot2::aes(y = upper, colour = model), linewidth = 0.12, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_line(data = spec, ggplot2::aes(colour = model), linewidth = 0.80, lineend = "round", na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.55, colour = "grey29") +
    ggplot2::facet_grid(response ~ impulse, scales = "free_y", switch = "y",
                        labeller = ggplot2::labeller(response = strip_labels)) +
    ggh4x::facetted_pos_scales(y = unname(y_scales[var_order])) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6), expand = ggplot2::expansion(mult = c(0.005, 0.01))) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Months", y = "") +
    ## base size 45 as in the published figures: it sets the strip margins and the padding between
    ## the outside strips and the axes; every visible text size is set explicitly below
    theme_publication(base_size = 45, base_family = base_family) +
    ggplot2::theme(
      panel.spacing     = grid::unit(2, "mm"),
      panel.spacing.x   = grid::unit(1, "mm"),
      panel.spacing.y   = grid::unit(4, "mm"),
      panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border      = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.title        = ggplot2::element_text(size = 55 * s, face = "plain"),
      axis.title.x      = ggplot2::element_text(angle = 0, vjust = 1, margin = ggplot2::margin(t = 2)),
      axis.title.y      = ggplot2::element_text(angle = 90, vjust = 1, margin = ggplot2::margin(r = -10, b = -5, l = 5)),
      axis.text         = ggplot2::element_text(size = 50 * s, colour = "black"),
      axis.text.y       = ggplot2::element_text(margin = ggplot2::margin(r = 2.5, l = -8, t = -5)),
      axis.text.x       = ggplot2::element_text(margin = ggplot2::margin(t = 1.2)),
      axis.line         = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.grid.major  = ggplot2::element_line(linewidth = 0.1, colour = "grey85"),
      panel.grid.minor  = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(colour = "white", fill = "white"),
      strip.text.y      = ggplot2::element_text(size = 55 * s, colour = "black", face = "plain", margin = ggplot2::margin(r = -1)),
      strip.text.x      = ggplot2::element_text(size = 50 * s, colour = "black", face = "plain"),
      strip.placement   = "outside",
      legend.position   = "none",
      plot.margin       = grid::unit(c(0, 0, 0, 0), "mm"),
      aspect.ratio      = 0.7
    )

  if (!is.null(labels)) {
    labels$response <- factor(labels$response, levels = var_order)
    labels$impulse  <- factor(impulse_label)
    p <- p + ggplot2::geom_text(data = labels, ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
                                size = 20 * s, fontface = "bold", colour = "black", family = base_family)
  }
  p
}

## ---- pairs-bootstrap bands (Appendix) -------------------------------------------------------

## The task1 objects store, for the point estimate and for every bootstrap draw, the VAR
## coefficient matrix B (regressors x equations: the constant first, then the p lag blocks of
## all variables) and the residual matrix. The climate variable is named "ACI" or "MOJO" in
## these matrices (or has no name); every column not in `macro_order` is taken as the climate index.

## ML covariance of the residuals (divides by n).
sigma_from_residuals <- function(U) stats::cov.wt(as.matrix(U), method = "ML")$cov

## Names of the equations in a saved matrix, with the climate variable called "climate".
var_names_from_matrix <- function(nms, macro_order) {
  nms <- if (is.null(nms)) rep("", length(macro_order) + 1L) else nms
  nms[is.na(nms) | !(nms %in% macro_order)] <- "climate"
  if (sum(nms == "climate") != 1L) stop("could not identify the climate variable among: ", paste(nms, collapse = ", "))
  nms
}

## Structural responses to the shock of the first variable in `var_order` (Cholesky
## identification) from a coefficient matrix B and a covariance Sigma, horizons 0..horizon.
## Returns an M x (horizon + 1) matrix with rows named by `var_order`; rows named in
## `cumulate` (logical named vector) are cumulated over the horizon.
irf_from_var_coefficients <- function(B, Sigma, p, var_order, horizon, cumulate, macro_order) {
  B <- as.matrix(B); M <- length(var_order)
  if (nrow(B) == M * p + 1L) B <- B[-1L, , drop = FALSE]        # drop the constant
  if (nrow(B) != M * p || ncol(B) != M) stop("B has dimension ", nrow(B), " x ", ncol(B), "; expected ", M * p, " x ", M)
  est_names <- var_names_from_matrix(colnames(B), macro_order)
  if (!is.null(rownames(B)) && !identical(var_names_from_matrix(rownames(B)[1:M], macro_order), est_names))
    warning("regressor order within a lag block differs from the equation order")
  perm <- match(var_order, est_names)                             # climate first, then the macro order
  A  <- lapply(seq_len(p), function(j) t(B[((j - 1L) * M + 1L):(j * M), , drop = FALSE])[perm, perm, drop = FALSE])
  Sig <- as.matrix(Sigma)[perm, perm, drop = FALSE]
  Fc <- matrix(0, M * p, M * p)                                   # companion matrix
  Fc[1:M, ] <- do.call(cbind, A)
  if (p > 1L) Fc[(M + 1L):(M * p), 1:(M * (p - 1L))] <- diag(M * (p - 1L))
  P <- t(chol(Sig))                                               # lower-triangular Cholesky factor
  ## J F^h J' P e_1 = J F^h (J' P e_1): propagate the companion state vector
  state <- c(P[, 1L], rep(0, M * (p - 1L)))
  out <- matrix(NA_real_, M, horizon + 1L, dimnames = list(var_order, as.character(0:horizon)))
  for (h in 0:horizon) {
    out[, h + 1L] <- state[1:M]
    state <- as.numeric(Fc %*% state)
  }
  for (v in var_order[cumulate[var_order]]) out[v, ] <- cumsum(out[v, ])
  out
}

## Point response and per-draw responses from a results/bootstrap/<id>_task1.rds object.
## Returns list(point = M x (H+1), draws = M x (H+1) x n_draws, n_draws, n_total).
pairs_bootstrap_irf <- function(task1, var_order, p, horizon, cumulate, macro_order) {
  irf_of <- function(B, U) irf_from_var_coefficients(B, sigma_from_residuals(U), p, var_order, horizon, cumulate, macro_order)
  point <- irf_of(task1$var_B, task1$var_resid)
  ok <- which(task1$boot$converged %in% TRUE &
              !vapply(task1$boot$var_B, is.null, logical(1)) & !vapply(task1$boot$var_resid, is.null, logical(1)))
  if (!length(ok)) stop("no converged bootstrap draw with saved VAR matrices in task1 object '", task1$id, "'")
  draws <- simplify2array(lapply(ok, function(i) irf_of(task1$boot$var_B[[i]], task1$boot$var_resid[[i]])))
  list(point = point, draws = draws, n_draws = length(ok), n_total = length(task1$boot$converged))
}

## Pointwise equal-tailed bands from the bootstrap draws. reverse_percentile = TRUE gives the
## basic (reverse-percentile) interval 2 * point - q_{1-a/2}, 2 * point - q_{a/2}; FALSE gives
## the percentile interval. Returns list(lower, upper), each M x (H+1).
pairs_bootstrap_bands <- function(point, draws, conf = 0.68, reverse_percentile = TRUE) {
  a <- 1 - conf
  q <- apply(draws, c(1, 2), stats::quantile, probs = c(a / 2, 1 - a / 2), names = FALSE, na.rm = TRUE)  # 2 x M x (H+1)
  if (reverse_percentile) {
    list(lower = 2 * point - q[2, , ], upper = 2 * point - q[1, , ])
  } else {
    list(lower = q[1, , ], upper = q[2, , ])
  }
}

## Copy of a band object whose lower/upper bands are replaced by the given M x (H+1) matrices
## (rows named by variable); the point response `irf` is kept.
replace_bands <- function(band_obj, lower, upper, source = "pairs bootstrap") {
  for (v in band_obj$var_order) {
    band_obj$bands[[v]]$lower <- as.numeric(lower[v, ])
    band_obj$bands[[v]]$upper <- as.numeric(upper[v, ])
  }
  band_obj$band_source <- source
  band_obj
}
