## Integrity and structure checks for the replication package: cleaned data, settings, the
## supplied result archives (manifest entries, dimensions and sizes), the point estimates, the
## figures and the summary-numbers file.
## run_all.R runs this after its last step; it can also be run alone from the package root:
##   Rscript code/funcs/validate_package.R
source("config.R")
source(file.path(ROOT, "code", "funcs", "load_data.R"))

failures <- character()
check_that <- function(ok, label) {
  if (length(ok) != 1L || is.na(ok) || !ok) failures <<- c(failures, label)
}
read_checked <- function(path, label) {
  if (!file.exists(path)) {
    failures <<- c(failures, paste0(label, " is missing"))
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) {
    failures <<- c(failures, paste0(label, " is unreadable: ", conditionMessage(e)))
    NULL
  })
}

required_data <- c(
  "aci_headline.csv", "aci_national_level.csv", "aci_national_transformed.csv",
  "aci_regional_level.csv", "aci_regional_transformed.csv", "labels_aci_national.csv",
  "labels_aci_regional.csv", "labels_noaa.csv", "macro_fredmd_2025-12.csv",
  "macro_transformations.csv", "nber_recessions.csv", "noaa_level.csv",
  "noaa_transformation_rules.csv", "noaa_transformed.csv", "README.md"
)
check_that(all(file.exists(file.path(PATHS$data, required_data))), "required cleaned data files")
check_that(identical(names(SPECS), c("comp_national", "comp_regional", "comp_noaa",
                                    "rank_national", "rank_regional", "rank_noaa")),
           "six specification identifiers")
check_that(identical(unname(vapply(SPECS, `[[`, integer(1), "grid_index")), c(8L, 4L, 12L, 1L, 6L, 27L)),
           "published CV grid indices")
check_that(CFG$start_date == "1975-01-01" && CFG$end_test_date == "2025-04-01" &&
             CFG$lags == 12L && CFG$irf_horizon == 36L && SEED == 123L,
           "declared sample, VAR, and seed settings")

expected_k <- c(comp_national = 5L, comp_regional = 34L, comp_noaa = 82L,
                rank_national = 5L, rank_regional = 34L, rank_noaa = 82L)
expected_converged <- c(comp_national = 999L, comp_regional = 1000L, comp_noaa = 1000L,
                        rank_national = 1000L, rank_regional = 1000L, rank_noaa = 1000L)

for (id in names(SPECS)) {
  cfg <- CFG
  cfg$climate_name <- "climate"
  d <- tryCatch(load_specification(SPECS[[id]], PATHS$data, cfg), error = function(e) {
    failures <<- c(failures, paste0("data-loader check for ", id, ": ", conditionMessage(e)))
    NULL
  })
  if (!is.null(d)) {
    check_that(nrow(d$comps) == 587L && ncol(d$comps) == expected_k[[id]], paste0(id, " panel dimensions"))
    check_that(length(d$dates_full) == 605L && length(d$dates) == 587L, paste0(id, " sample dates"))
  }

  pe <- read_checked(file.path(PATHS$estimates, paste0(id, ".rds")), paste0(id, " point estimate"))
  if (!is.null(pe)) {
    check_that(identical(pe$id, id) && length(pe$weights) == expected_k[[id]], paste0(id, " point-estimate shape"))
    check_that(pe$lags == 12L && pe$end_test_date == "2025-04-01", paste0(id, " point-estimate settings"))
    lambda_grid <- EST$lambda_grid(SPECS[[id]]$rank)
    check_that(any(abs(as.numeric(pe$lambda) - lambda_grid) <= 1e-10 * pmax(1, abs(lambda_grid))),
               paste0(id, " point-estimate lambda is on the CV grid"))
    if (is.null(pe$cv))
      check_that(isTRUE(all.equal(as.numeric(pe$lambda), SPECS[[id]]$lambda, tolerance = 1e-12)),
                 paste0(id, " supplied point-estimate lambda"))
  }

  t1 <- read_checked(file.path(PATHS$bootstrap, paste0(id, "_task1.rds")), paste0(id, " task 1"))
  if (!is.null(t1)) {
    check_that(identical(t1$id, id) && t1$origin == 216L, paste0(id, " task-1 identity"))
    check_that(identical(dim(t1$boot$weights), c(1000L, expected_k[[id]])), paste0(id, " task-1 dimensions"))
    check_that(t1$n_draws == 1000L && t1$n_converged == expected_converged[[id]], paste0(id, " task-1 draw counts"))
    check_that(t1$sample_n == 560L && t1$block_length == 13L, paste0(id, " task-1 sample settings"))
    check_that(isTRUE(all.equal(as.numeric(t1$lambda), SPECS[[id]]$lambda, tolerance = 5e-7)), paste0(id, " task-1 lambda"))
    check_that(is.null(t1$session$seed) || t1$session$seed == 123L, paste0(id, " task-1 seed"))
  }

  t2 <- read_checked(file.path(PATHS$bootstrap, paste0(id, "_task2.rds")), paste0(id, " task 2"))
  if (!is.null(t2)) {
    check_that(identical(t2$id, id), paste0(id, " task-2 identity"))
    check_that(identical(dim(t2$W_mat), c(216L, expected_k[[id]])), paste0(id, " task-2 dimensions"))
    check_that(identical(as.integer(t2$origins), 1:216), paste0(id, " task-2 origins"))
    check_that(isTRUE(all.equal(as.numeric(t2$lambda), SPECS[[id]]$lambda, tolerance = 5e-7)), paste0(id, " task-2 lambda"))
    check_that(is.null(t2$session$seed) || t2$session$seed == 123L, paste0(id, " task-2 seed"))
  }
}

for (id in c(names(SPECS), "benchmark")) {
  irf <- read_checked(file.path(PATHS$irf_bands, paste0(id, ".rds")), paste0(id, " IRF bands"))
  if (!is.null(irf)) {
    check_that(identical(dim(irf$irf_point), c(37L, 7L)) && length(irf$bands) == 7L,
               paste0(id, " IRF dimensions"))
    check_that(irf$n_obs == 585L && irf$block_length == 8L, paste0(id, " IRF sample settings"))
    check_that(irf$B_outer == 1000L && irf$B_inner == 100L && irf$B_se == 200L && irf$seed == 123L,
               paste0(id, " IRF bootstrap settings"))
    check_that(all(vapply(irf$bands, function(z) all(c("irf", "lower", "upper", "cumulated") %in% names(z)) &&
                                      length(z$irf) == 37L && length(z$lower) == 37L && length(z$upper) == 37L,
                              logical(1))), paste0(id, " IRF band fields"))
  }
}

manifest_file <- file.path(PATHS$bootstrap, "manifest.csv")
if (file.exists(manifest_file)) {
  manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
  check_that(nrow(manifest) == 12L && identical(sort(manifest$filename),
             sort(paste0(rep(names(SPECS), each = 2L), "_task", 1:2, ".rds"))),
             "bootstrap manifest entries")
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(PATHS$bootstrap, manifest$filename[i])
    if (file.exists(path))
      check_that(file.info(path)$size == manifest$size_bytes[i], paste0(manifest$filename[i], " size"))
  }
} else failures <- c(failures, "bootstrap manifest is missing")

## the 38 figure files of the package (list of figures of the paper): 5 per-specification families x 6
## specifications, 6 main-text series and weight figures, the NOAA rank-usage figure, and the
## static legend strip
expected_figures <- c(
  paste0(rep(c("irf_", "irfboot_", "boot_weights_", "boot_index_", "expanding_weights_"),
             each = length(SPECS)), rep(names(SPECS), times = 5L), ".png"),
  paste0(c("fitted_aci", "fitted_noaa", "weights_comp_aci", "weights_rank_aci",
           "weights_comp_noaa", "weights_rank_noaa", "rank_usage_noaa"), ".png"),
  file.path("static", "Naming_ACI.png")
)
check_that(length(expected_figures) == 38L && !anyDuplicated(expected_figures), "38 unique expected figures")
figure_paths <- file.path(PATHS$figures, expected_figures)
figure_ok <- file.exists(figure_paths) & !is.na(file.info(figure_paths)$size) & file.info(figure_paths)$size > 0
check_that(all(figure_ok), paste0("figure files missing or empty: ",
                                  paste(expected_figures[!figure_ok], collapse = ", ")))

## the in-text numbers written by step 6
summary_file <- file.path(PATHS$summary, "summary_numbers.txt")
check_that(file.exists(summary_file) && isTRUE(file.info(summary_file)$size > 0),
           "results/summary_numbers.txt missing or empty")

if (length(failures)) {
  cat("PACKAGE VALIDATION FAILED\n", paste0("- ", unique(failures), collapse = "\n"), "\n", sep = "")
  quit(status = 1L)
}
cat("PACKAGE VALIDATION PASSED: data, settings, 6 point estimates, 7 IRF band archives, 12 stability archives, 38 figures, and the summary numbers\n")
