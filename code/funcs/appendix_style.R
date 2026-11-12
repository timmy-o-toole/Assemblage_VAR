## Display order and named colors of the appendix weight and NOAA-usage figures.
## These mappings affect presentation only; numerical summaries are prepared elsewhere.

APPENDIX_ACI_ORDER <- c("WP90", "Sea.Level", "Rx5Day", "T10", "T90")
APPENDIX_NOAA_BOOT_ORDER <- c("Heat.DD", "Palmer.Z", "Precip.", "Temp.", "Others")
APPENDIX_NOAA_EXPANDING_ORDER <- c("Temp.", "Precip.", "Palmer.Z", "Heat.DD", "Others")
APPENDIX_RANK_ORDER <- paste0("Q", 1:5)
APPENDIX_REGION_ORDER <- c("Up.MidW.", "Ohio.V.", "N.West", "N.Rockies", "N.East", "Rest")

APPENDIX_FAMILY_COLOURS <- c(
  WP90 = "#001A43", Sea.Level = "#328CC1", Rx5Day = "#1D7874",
  T10 = "#F4A261", T90 = "#8C1D40", Heat.DD = "#F4A261",
  Palmer.Z = "#1D7874", Precip. = "#328CC1", Temp. = "#8C1D40",
  Others = "#6B705C"
)
APPENDIX_RANK_COLOURS <- c(
  Q1 = "#001330", Q2 = "#1D3F6E", Q3 = "#3C68BB",
  Q4 = "#7FA8D9", Q5 = "#C9DCF0"
)
APPENDIX_REGION_COLOURS <- c(
  "Up.MidW." = "#D44233", "Ohio.V." = "#197A81", "N.West" = "#FF8500",
  "N.Rockies" = "#3C68BB", "N.East" = "#EFAE0E", "Rest" = "#782F26"
)
APPENDIX_EXPANDING_ALPHA <- 0.88
APPENDIX_USAGE_ALPHA <- 0.90

appendix_component_levels <- function(id, figure = c("bootstrap", "expanding")) {
  figure <- match.arg(figure)
  if (id %in% c("comp_national", "comp_regional")) return(APPENDIX_ACI_ORDER)
  if (id == "comp_noaa") {
    if (figure == "bootstrap") return(APPENDIX_NOAA_BOOT_ORDER)
    return(APPENDIX_NOAA_EXPANDING_ORDER)
  }
  stop("no component display order for ", id, call. = FALSE)
}

appendix_display_levels <- function(present, order) {
  present <- unique(as.character(present))
  if (!setequal(present, order))
    stop("display groups do not match the expected paper categories: ",
         paste(setdiff(present, order), collapse = ", "), call. = FALSE)
  order
}

appendix_usage_levels <- function(present_keys) {
  keys <- unlist(lapply(APPENDIX_NOAA_BOOT_ORDER, function(family)
    paste(family, APPENDIX_REGION_ORDER, sep = "::")), use.names = FALSE)
  if (!all(present_keys %in% keys)) stop("unexpected NOAA region-family cell", call. = FALSE)
  keys[keys %in% present_keys]
}
