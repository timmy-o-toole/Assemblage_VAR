## =============================================================================
## code/funcs/theme.R -- shared ggplot2 theme and colours for all figures.
## Built on ggthemes::theme_foundation.
## Fonts: "STIX Two Text" if it is installed locally, otherwise the package
## falls back to the default serif family, so no network access is needed.
## =============================================================================

figure_font <- function() {
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character())
    for (f in c("STIX Two Text", "STIXGeneral", "Times New Roman")) if (f %in% fams) return(f)
  }
  "serif"
}

theme_publication <- function(base_size = 30, base_family = figure_font()) {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3), hjust = 0.5),
      text = ggplot2::element_text(),
      panel.background = ggplot2::element_rect(colour = NA),
      plot.background = ggplot2::element_rect(colour = NA),
      panel.border = ggplot2::element_rect(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(angle = 45, vjust = -0.2),
      axis.text = ggplot2::element_text(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.ticks = ggplot2::element_line(),
      panel.grid.major = ggplot2::element_line(colour = "#f0f0f0"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.6, "cm"),
      legend.margin = ggplot2::margin(-20, 5, 5, -20),
      legend.box.margin = ggplot2::margin(-5, -5, -5, -5),
      legend.title = ggplot2::element_text(face = "italic"),
      plot.margin = grid::unit(c(5, 5, 5, 5), "mm"),
      strip.background = ggplot2::element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

## Human-readable names used in figure titles and legends
SPEC_TITLES <- c(comp_national = "National", comp_regional = "Regional", comp_noaa = "NOAA",
                 rank_national = "National", rank_regional = "Regional", rank_noaa = "NOAA")
VAR_LABELS  <- c(climate = "Climate index", INDPRO = "Industrial production", PCEPI = "PCE inflation",
                 UNRATE = "Unemployment rate", SHADOWRATE = "Shadow rate", EXUSUKx = "USD/GBP exchange rate",
                 HOUST = "Housing starts")

## Save helper: PNG via ragg if available (crisper text), else grDevices::png through ggsave.
save_figure <- function(plot, file, width, height, units = "mm", dpi = 300, scale = 1) {
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(file, plot = plot, width = width, height = height, units = units, dpi = dpi, scale = scale, device = dev, bg = "white")
  invisible(file)
}
