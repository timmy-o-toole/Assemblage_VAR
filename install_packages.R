## install_packages.R -- install the CRAN packages needed by the replication package and check
## their versions. The tested versions are listed in requirements.csv. Newer CRAN versions are
## expected to work; older ones are rejected (for example, the figures need ggplot2 >= 3.5.0).
##   Rscript install_packages.R        (run once, from the package root; needs internet access
##                                      only if a package is missing)
req <- read.csv("requirements.csv", stringsAsFactors = FALSE)
required <- req$package[!grepl("optional", req$role, ignore.case = TRUE)]
optional <- setdiff(req$package, required)
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) {
  cat("installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}
installed <- rownames(installed.packages())
still_missing <- setdiff(required, installed)
if (length(still_missing))
  stop("required packages could not be installed: ", paste(still_missing, collapse = ", "), call. = FALSE)
too_old <- character()
for (p in required) {
  tested <- req$version[req$package == p]
  if (packageVersion(p) < package_version(tested))
    too_old <- c(too_old, sprintf("%s %s (tested with %s)", p, packageVersion(p), tested))
}
if (length(too_old))
  stop("packages older than the tested versions, please update them:\n  ",
       paste(too_old, collapse = "\n  "), call. = FALSE)
cat("all required packages are installed at or above the tested versions\n")
optional_missing <- setdiff(optional, installed)
## the role in requirements.csv says what each missing optional package affects
if (length(optional_missing))
  cat("optional packages not installed:\n",
      sprintf("  %-12s %s\n", optional_missing, req$role[match(optional_missing, req$package)]), sep = "")
