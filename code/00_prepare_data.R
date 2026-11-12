## =============================================================================
## 00_prepare_data.R -- Step 0 (optional): rebuild the cleaned data files in data/
## from the raw source files in data/raw/.
##
## Not called by run_all.R: the package ships the cleaned files, and this step
## documents how they were built. It reads the raw files, seasonally adjusts the
## climate series with X-13ARIMA-SEATS (RJDemetra, needs Java), selects and
## transforms the series, writes every file to data/_prepared/ (the shipped files
## in data/ are never overwritten) and compares each rebuilt file with the
## shipped one. The transformed panels are built with the loader's own functions
## (code/funcs/load_data.R). Raw files and where to obtain them: data/raw/README.md (sizes and
## SHA-256 checksums: data/raw/MANIFEST.csv).
##
## Usage:  Rscript code/00_prepare_data.R    (from the package root; 5-10 minutes)
## =============================================================================
source("config.R")
source(file.path(ROOT, "code", "funcs", "load_data.R"))   # transformations shared with the loader

for (p in c("readxl", "RJDemetra", "zoo"))
  if (!requireNamespace(p, quietly = TRUE))
    stop("package '", p, "' is needed for the data preparation (see requirements.csv; ",
         "RJDemetra needs rJava and a Java 8+ runtime)", call. = FALSE)

RAW <- file.path(PATHS$data, "raw")
OUT <- file.path(PATHS$data, "_prepared")
TOL <- 1e-8                                     # maximal relative difference to the shipped files
ACI_WINDOW    <- c("1961-01-01", "2025-05-01")  # months used from the ACI workbook
NOAA_WINDOW   <- c("1961-01-01", "2025-12-01")  # months used from the NOAA export
SHADOW_WINDOW <- c("1990-04-01", "2022-02-01")  # months in which SHADOWRATE is the Wu-Xia shadow rate
US_REGIONS    <- c("ALA", "CEA", "CWP", "MID", "SEA", "SPL", "SWP")   # the seven U.S. ACI regions

raw_files <- c(
  "ACI_Data_Site.xlsx"         = "Actuaries Climate Index workbook",
  "noaa_macrobond_export.xlsx" = "Macrobond export of NOAA regional climate series and Mauna Loa CO2",
  "2025-12-md.csv"             = "FRED-MD vintage 2025-12",
  "WuXiaShadowRate.xlsx"       = "Wu-Xia shadow federal funds rate")
raw_missing <- !file.exists(file.path(RAW, names(raw_files)))
if (any(raw_missing))
  stop("raw file(s) not found in data/raw/:\n",
       paste0("  - ", names(raw_files)[raw_missing], " (", raw_files[raw_missing], ")", collapse = "\n"),
       "\nSee data/raw/README.md for their sources and how to obtain them.", call. = FALSE)
f_aci   <- file.path(RAW, "ACI_Data_Site.xlsx")
f_noaa  <- file.path(RAW, "noaa_macrobond_export.xlsx")
f_fred  <- file.path(RAW, "2025-12-md.csv")
f_wuxia <- file.path(RAW, "WuXiaShadowRate.xlsx")
dir.create(OUT, showWarnings = FALSE)
t_start <- Sys.time()

## ---- helpers ------------------------------------------------------------------------------

## CSV with LF line endings on every platform (binary connection), as the shipped files
write_lf <- function(df, fname, ...) {
  con <- file(file.path(OUT, fname), "wb"); on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE, ...)
}
write_panel <- function(x, name) {
  df <- cbind(date = rownames(x), as.data.frame(x, check.names = FALSE))
  write_lf(df, name, na = "")
  cat(sprintf("  wrote %-32s %4d x %d\n", name, nrow(df), ncol(df) - 1L))
}
write_table <- function(df, name) {
  write_lf(df, name)
  cat(sprintf("  wrote %-32s %4d x %d\n", name, nrow(df), ncol(df)))
}
in_window <- function(dates, window) dates >= window[1] & dates <= window[2]
check_monthly <- function(dates, what) {
  if (!identical(dates, format(seq(as.Date(dates[1]), by = "month", length.out = length(dates)))))
    stop(what, ": dates are not a gap-free monthly sequence", call. = FALSE)
}

## X-13ARIMA-SEATS with the RJDemetra default specification (RSA5c); returns the SA series.
x13_sa <- function(x, start) {
  fit <- RJDemetra::x13(stats::ts(as.numeric(x), start = start, frequency = 12))
  as.numeric(fit$final$series[, "sa"])
}
x13_panel <- function(m, label) {
  start <- as.integer(strsplit(rownames(m)[1], "-")[[1]][1:2])
  out <- m
  for (j in seq_len(ncol(m))) out[, j] <- x13_sa(m[, j], start)
  cat(sprintf("  X-13: %s, %d series\n", label, ncol(m)))
  out
}

## ---- ACI: national and regional components ---------------------------------------------------

## One "<component> Unstandardized" sheet of the ACI workbook -> months x regions. The block
## "Unstandardized: Monthly Unsmoothed" has one row per region; years and months are in the
## third and second row above its title.
read_aci_sheet <- function(sheet) {
  x <- as.data.frame(suppressMessages(readxl::read_excel(f_aci, sheet = sheet, col_names = FALSE,
                                                          .name_repair = "minimal")))
  lab <- x[[1]]
  i0 <- which(lab == "Unstandardized: Monthly Unsmoothed")
  if (length(i0) != 1) stop("ACI sheet '", sheet, "': monthly unsmoothed block not found", call. = FALSE)
  yr <- suppressWarnings(as.integer(unlist(x[i0 - 3, -1])))
  mo <- suppressWarnings(as.integer(unlist(x[i0 - 2, -1])))
  ok <- !is.na(yr) & !is.na(mo)
  blank <- which(is.na(lab) & seq_along(lab) > i0)
  rows <- (i0 + 1):(if (length(blank)) blank[1] - 1 else length(lab))
  rows <- rows[grepl("^[A-Z]{3}$", lab[rows])]                 # skips "No data for ... region"
  m <- vapply(rows, function(r) as.numeric(unlist(x[r, -1]))[ok], numeric(sum(ok)))
  dimnames(m) <- list(sprintf("%04d-%02d-01", yr[ok], mo[ok]), lab[rows])
  check_monthly(rownames(m), paste("ACI sheet", sheet))
  m[in_window(rownames(m), ACI_WINDOW), , drop = FALSE]
}

cat("ACI workbook\n")
aci_sheets <- c(Sea.Level = "Sea Level Unstandardized", CDD = "CDD Unstandardized",
                Rx5Day = "Rx5day Unstandardized", T10 = "T10 Unstandardized",
                T90 = "T90 Unstandardized", WP90 = "WP90 Unstandardized")
aci <- lapply(aci_sheets, read_aci_sheet)
if (!all(vapply(aci, function(m) identical(rownames(m), rownames(aci[[1]])), logical(1))))
  stop("ACI workbook: the component sheets cover different months", call. = FALSE)
## national: the USA row of all six components (CDD is kept for the level file and dropped later)
aci_nat <- sapply(aci, function(m) m[, "USA"])
## regional: the seven U.S. regions of the five components without CDD (no sea level for MID)
aci_reg <- do.call(cbind, lapply(setdiff(names(aci), "CDD"), function(cp) {
  m <- aci[[cp]][, intersect(US_REGIONS, colnames(aci[[cp]])), drop = FALSE]
  colnames(m) <- paste0(cp, ".", colnames(m)); m
}))
write_panel(x13_panel(aci_nat, "ACI national"), "aci_national_level.csv")
write_panel(x13_panel(aci_reg, "ACI regional"), "aci_regional_level.csv")

## ---- NOAA panel ---------------------------------------------------------------------------

cat("NOAA export\n")
x <- as.data.frame(suppressMessages(readxl::read_excel(f_noaa, sheet = 1, col_names = FALSE,
                                                        .name_repair = "minimal")))
nm <- as.character(unlist(x[1, -1]))          # first row: Macrobond series descriptions
x <- x[-1, , drop = FALSE]
## dates: "dd.mm.yyyy" text before 1900, Excel serial numbers afterwards
d <- x[[1]]
dates <- ifelse(grepl("^\\d{2}\\.\\d{2}\\.\\d{4}$", d),
                format(as.Date(d, format = "%d.%m.%Y")),
                format(as.Date(suppressWarnings(as.numeric(d)), origin = "1899-12-30")))
check_monthly(dates, "NOAA export")
if (anyNA(nm)) stop("NOAA export: missing series name in the first row", call. = FALSE)
keep <- !duplicated(nm)                        # repeated series: first occurrence kept
raw <- sapply(x[-1][keep], as.numeric)
dimnames(raw) <- list(dates, nm[keep])
raw <- raw[in_window(dates, NOAA_WINDOW), , drop = FALSE]
cat(sprintf("  %d series (%d repeated columns dropped), %s to %s\n",
            ncol(raw), sum(!keep), rownames(raw)[1], rownames(raw)[nrow(raw)]))
raw[] <- apply(raw, 2, zoo::na.approx, maxgap = 2, na.rm = FALSE)   # gaps of at most two months

## X-13 where the series allows it: at least five years, no missing values, non-zero variance and
## a successful fit; otherwise the series is left unadjusted.
start <- as.integer(strsplit(NOAA_WINDOW[1], "-")[[1]][1:2])
adjusted <- logical(ncol(raw))
sa <- raw
for (j in seq_len(ncol(raw))) {
  v <- raw[, j]
  if (length(v) < 60 || anyNA(v) || stats::sd(v) == 0) next
  s <- tryCatch(x13_sa(v, start), error = function(e) NULL)
  if (!is.null(s)) { sa[, j] <- s; adjusted[j] <- TRUE }
}
cat(sprintf("  X-13: NOAA, %d of %d series adjusted, %d left unadjusted\n",
            sum(adjusted), ncol(raw), sum(!adjusted)))

## Transformation rule per series: log difference for strictly positive series, first difference
## otherwise ("zero" for series with fewer than ten observations or no variation; these are dropped
## below).
noaa_rule <- function(v) {
  z <- stats::na.omit(v)
  if (length(z) < 10 || stats::sd(z) == 0) "zero" else if (all(z > 0)) "logdiff" else "diff"
}
rules <- vapply(colnames(sa), function(cn) noaa_rule(sa[, cn]), character(1))

## Series selection: transform, drop the first month, set non-finite values to missing, drop series
## with more than 7.5% missing values or zero variance, interpolate internal gaps, and keep the
## series without remaining (leading or trailing) gaps.
st <- vapply(colnames(sa), function(cn) transform_noaa_series(sa[, cn], rules[[cn]]), numeric(nrow(sa)))
st <- st[-1, , drop = FALSE]
st[!is.finite(st)] <- NA
sdv <- apply(st, 2, stats::sd, na.rm = TRUE)
ok <- colMeans(is.na(st)) <= 0.075 & !is.na(sdv) & sdv > 0
ok[ok] <- apply(st[, ok, drop = FALSE], 2, function(v) !anyNA(zoo::na.approx(v, na.rm = FALSE)))
cat(sprintf("  %d of %d series kept (%s; %d seasonally adjusted)\n", sum(ok), length(ok),
            paste(names(table(rules[ok])), table(rules[ok]), collapse = ", "), sum(adjusted[ok])))
noaa_cols <- colnames(sa)[ok]
write_panel(sa[, noaa_cols, drop = FALSE], "noaa_level.csv")
write_table(data.frame(name = noaa_cols, rule = unname(rules[noaa_cols])), "noaa_transformation_rules.csv")

## ---- transformed panels (the loader's functions, applied to the rebuilt level files) ----------

cat("Transformed panels\n")
nat <- build_aci_national(OUT, check = FALSE)
write_panel(nat$comps, "aci_national_transformed.csv")
write_panel(data.frame(ACI = nat$headline, row.names = names(nat$headline)), "aci_headline.csv")
write_panel(build_aci_regional(OUT, check = FALSE), "aci_regional_transformed.csv")
write_panel(build_noaa(OUT, check = FALSE), "noaa_transformed.csv")

## ---- labels ---------------------------------------------------------------------------------

cat("Labels\n")
write_table(data.frame(index = seq_along(colnames(nat$comps)), name = colnames(nat$comps)), "labels_aci_national.csv")
reg_cols <- colnames(aci_reg)
write_table(data.frame(index = seq_along(reg_cols), name = reg_cols, variable = sub("\\.[A-Z]{3}$", "", reg_cols),
                       region = sub("^.*\\.([A-Z]{3})$", "\\1", reg_cols)), "labels_aci_regional.csv")
## NOAA names: "<country>, <region>, Weather & Atmosphere, <variable>" or "<region>, <variable>".
## Aggregated groups as in Tables C.1 and C.2 of the paper; all other regions are "Rest" and all
## other variables "Others".
reg_map <- c("Upper Midwest" = "Up.MidW.", "Ohio Valley" = "Ohio.V.", "Northwest" = "N.West",
             "Northern Rockies & Plains" = "N.Rockies", "Northeast" = "N.East")
var_map <- c("Temperature" = "Temp.", "Temperature, Average" = "Temp.", "Precipitation" = "Precip.",
             "Heating Degree Days" = "Heat.DD", "Palmer Z-Index, Index" = "Palmer.Z")
parts <- strsplit(noaa_cols, ", ", fixed = TRUE)
wa <- vapply(parts, function(p) length(p) >= 3 && p[3] == "Weather & Atmosphere", logical(1))
region_raw   <- ifelse(wa, vapply(parts, `[`, "", 2), vapply(parts, `[`, "", 1))
variable_raw <- vapply(seq_along(parts), function(i)
  paste(parts[[i]][-seq_len(if (wa[i]) 3 else 1)], collapse = ", "), character(1))
write_table(data.frame(index = seq_along(noaa_cols), name = noaa_cols, region_raw = region_raw, variable_raw = variable_raw,
                       region_agg = ifelse(region_raw %in% names(reg_map), reg_map[region_raw], "Rest"),
                       variable_agg = ifelse(variable_raw %in% names(var_map), var_map[variable_raw], "Others")),
            "labels_noaa.csv")

## ---- macro block: FRED-MD and the Wu-Xia shadow rate ---------------------------------------------

cat("FRED-MD\n")
macro_series <- c("INDPRO", "PCEPI", "UNRATE", "SHADOWRATE", "EXUSUKx", "HOUST")
fred <- utils::read.csv(f_fred, stringsAsFactors = FALSE)
tcodes <- fred[1, ]                               # row 1: FRED-MD transformation codes
fred <- fred[2:(nrow(fred) - 1), ]                # final month (2025-11, incomplete) dropped
rownames(fred) <- format(as.Date(fred$sasdate, format = "%m/%d/%Y"))
check_monthly(rownames(fred), "FRED-MD")
## SHADOWRATE: the Wu-Xia shadow federal funds rate (sheet "Data" of the Atlanta Fed workbook) in
## SHADOW_WINDOW, the FRED-MD effective federal funds rate FEDFUNDS before and after (its
## transformation code is that of FEDFUNDS)
wx <- as.data.frame(suppressMessages(readxl::read_excel(f_wuxia, sheet = "Data")))
wx_col <- grep("^Wu-Xia shadow federal funds rate", names(wx))
if (length(wx_col) != 1) stop("Wu-Xia workbook: shadow rate column not found in sheet 'Data'", call. = FALSE)
win <- in_window(rownames(fred), SHADOW_WINDOW)
wx_rate <- wx[[wx_col]][match(rownames(fred)[win], format(as.Date(wx[[1]])))]
if (anyNA(wx_rate)) stop("Wu-Xia workbook: shadow rate missing for months in ",
                         paste(SHADOW_WINDOW, collapse = " to "), call. = FALSE)
fred$SHADOWRATE <- replace(fred$FEDFUNDS, win, wx_rate)
tcodes$SHADOWRATE <- tcodes$FEDFUNDS
cat(sprintf("  SHADOWRATE: Wu-Xia shadow rate %s to %s (%d months), FEDFUNDS otherwise\n",
            substr(SHADOW_WINDOW[1], 1, 7), substr(SHADOW_WINDOW[2], 1, 7), sum(win)))
write_panel(fred[, macro_series], "macro_fredmd_2025-12.csv")
write_table(data.frame(series = macro_series, fredmd_tcode = as.integer(unlist(tcodes[macro_series])),
                       transformation_used = c("100*dlog", "100*dlog", "diff", "level", "100*dlog", "100*dlog")),
            "macro_transformations.csv")

## ---- NBER recessions ----------------------------------------------------------------------------

## NBER business cycle reference dates (peak and trough months), nine recessions 1960-2020, from
## https://www.nber.org/research/data/us-business-cycle-expansions-and-contractions
cat("NBER recessions\n")
nber <- data.frame(
  start = as.Date(c("1960-04-01", "1969-12-01", "1973-11-01", "1980-01-01", "1981-07-01",
                    "1990-07-01", "2001-03-01", "2007-12-01", "2020-02-01")),
  end   = as.Date(c("1961-02-01", "1970-11-01", "1975-03-01", "1980-07-01", "1982-11-01",
                    "1991-03-01", "2001-11-01", "2009-06-01", "2020-04-01")))
write_table(nber, "nber_recessions.csv")

## ---- compare with the shipped files ----------------------------------------------------------

## X-13ARIMA-SEATS (JDemetra+, run in Java) gives the same result on repeated runs on one machine,
## but its automatic model selection and estimation can differ slightly across platforms and Java
## versions. Differences in the files that depend on it are therefore reported, not treated as
## errors; every other file must match the shipped one.
x13_files <- c("aci_national_level.csv", "aci_regional_level.csv", "noaa_level.csv", "aci_national_transformed.csv",
               "aci_headline.csv", "aci_regional_transformed.csv", "noaa_transformed.csv")
other_files <- c("noaa_transformation_rules.csv", "labels_aci_national.csv", "labels_aci_regional.csv", "labels_noaa.csv",
                 "macro_fredmd_2025-12.csv", "macro_transformations.csv", "nber_recessions.csv")
compare_file <- function(name) {
  a <- utils::read.csv(file.path(OUT, name), check.names = FALSE, stringsAsFactors = FALSE)
  b <- utils::read.csv(file.path(PATHS$data, name), check.names = FALSE, stringsAsFactors = FALSE)
  num <- vapply(b, is.numeric, logical(1))
  res <- data.frame(file = name, identical = unname(tools::md5sum(file.path(OUT, name)) == tools::md5sum(file.path(PATHS$data, name))),
                    same_shape = FALSE, cols_above_tol = NA_integer_, n_cols = sum(num), max_abs = NA_real_, max_rel = NA_real_)
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) return(res)
  res$same_shape <- identical(a[!num], b[!num]) && identical(is.na(as.matrix(a[num])), is.na(as.matrix(b[num])))
  dd <- abs(as.matrix(a[num]) - as.matrix(b[num])); rel <- dd / pmax(1, abs(as.matrix(b[num])))
  res$max_abs <- max(c(0, dd), na.rm = TRUE); res$max_rel <- max(c(0, rel), na.rm = TRUE)
  res$cols_above_tol <- sum(apply(rel, 2, max, na.rm = TRUE) > TOL)
  res
}
cmp <- do.call(rbind, lapply(c(x13_files, other_files), compare_file))
cmp$status <- ifelse(cmp$identical, "identical",
              ifelse(cmp$same_shape & cmp$max_rel <= TOL, "match",
              ifelse(cmp$same_shape & cmp$file %in% x13_files, "differs (X-13)", "MISMATCH")))
cat(sprintf("\nComparison with the shipped files in data/ (tolerance %g on |rebuilt - shipped| / max(1, |shipped|))\n", TOL))
cat(sprintf("  %-30s %-15s %11s %10s %10s\n", "file", "status", "cols > tol", "max abs", "max rel"))
cat(sprintf("  %-30s %-15s %5s of %3d %10.2e %10.2e\n", cmp$file, cmp$status,
            ifelse(is.na(cmp$cols_above_tol), "-", cmp$cols_above_tol), cmp$n_cols, cmp$max_abs, cmp$max_rel), sep = "")
cat(sprintf("run time %.1f min, rebuilt files in data/_prepared/\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
if (any(cmp$status == "MISMATCH"))
  stop("rebuilt files do not match the shipped data: ", paste(cmp$file[cmp$status == "MISMATCH"], collapse = ", "), call. = FALSE)
if (any(cmp$status == "differs (X-13)"))
  cat("note: the seasonally adjusted series differ from the shipped ones beyond the tolerance for the",
      "columns counted above (X-13 results depend on the platform); all other files match.\n")
