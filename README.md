# Building Macroeconomically Relevant Climate Indices via the Assemblage VAR

**Christophe Barrette, Philippe Goulet Coulombe and Tim Reinicke**
*The Econometrics Journal*, forthcoming.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23976450.svg)](https://doi.org/10.5281/zenodo.23976450)
![R >= 4.3](https://img.shields.io/badge/R-%3E%3D%204.3-blue)
[![License: MIT](https://img.shields.io/badge/code-MIT-green.svg)](LICENSE)

[Paper](https://arxiv.org/abs/2610.10043) · [Replication package on Zenodo](https://doi.org/10.5281/zenodo.23976450)

The Assemblage VAR estimates nonnegative aggregation weights on the subcomponents of a composite
index jointly with the VAR, so that the assembled index is the one most relevant for the
macroeconomic variables. It comes in two variants, component space and rank space. This repository
contains the code and cleaned data to rerun the paper's estimates, vary them, and apply the method
to other data.

## Quick start

From the repository root, with R 4.3 or newer:

```bash
Rscript install_packages.R   # once
Rscript run_all.R            # ~5-10 min; all estimates and figures of the paper
```

On Windows without `Rscript` on the PATH, `source()` the two files from R instead.

## Vary the settings

All settings are in [`config.R`](config.R). Rerun one specification with
`Rscript code/01_estimate_weights.R comp_national`. The result,
`results/point_estimates/comp_national.rds`, contains `$weights` and `$index`.

| To change | In `config.R` |
|---|---|
| Penalty | `SPECS$<spec>$lambda`, or `RUN_CV=TRUE` to select it by cross-validation |
| Sample, lags, ordering, horizon | `CFG$start_date`, `CFG$end_test_date`, `CFG$lags`, `CFG$macro_order`, `CFG$irf_horizon` |
| Panel and variant | the six `SPECS`: ACI national (5 subcomponents), ACI regional (34), NOAA (82), each in component and rank space |

Impulse responses for changed settings, with a quick bootstrap:
`B_OUTER=20 B_INNER=5 B_SE=5 RUN_IRF_BOOTSTRAP=TRUE Rscript code/02_irf_bands.R comp_national`.

## Use your own data

You need two monthly CSV files, each with a `date` column (`YYYY-MM-01`) followed by one column per
series:

- **Subcomponents** of the composite index, one column each, transformed to stationarity and
  standardised (as in `data/aci_national_transformed.csv`).
- **Macroeconomic variables** of the VAR, transformed to stationarity, in Cholesky order. The
  assembled index enters first.

```r
source("config.R"); source("code/funcs/load_data.R"); source("code/funcs/estimation.R"); source("code/funcs/cv.R")

comps <- read_panel("my_subcomponents.csv")
macro <- read_panel("my_macro.csv")
dates <- intersect(rownames(comps), rownames(macro))
comps <- comps[dates, , drop = FALSE]                    # rank space: comps <- rank_sort(comps)$sorted
d <- list(comps = comps, macro = macro[dates, ], dates = dates, end_test_date = tail(dates, 1),
          climate_name = "climate",
          bench = data.frame(climate = rowMeans(comps), macro[dates, ], check.names = FALSE))  # fixed-weight benchmark

spec <- list(rank = FALSE, level = 1L)
est  <- prepare_estimation(d, p = 12, holdout = 0)
lam  <- cross_validate(est, d, spec, EST$lambda_grid(spec$rank), EST)$lambda
fit  <- solve_weights(est, d, spec, lam, EST)

round(sort(fit$weights, decreasing = TRUE), 3)          # aggregation weights
index <- assembled_index(d, fit$weights)                 # assembled index
```

The assembled index can then replace the fixed-weight benchmark in any VAR
(see `code/02_irf_bands.R` for the impulse responses used in the paper).

## Replication and data

For exact replication, use the archived package on
[Zenodo](https://doi.org/10.5281/zenodo.23976450), which also contains the raw data and full
documentation. Data sources are described in [`data/README.md`](data/README.md) and
[`data/raw/README.md`](data/raw/README.md). The MIT licence covers the code; the data remain under
their providers' terms.

## Citation

```bibtex
@article{BarretteGouletCoulombeReinicke2026,
  title   = {Building Macroeconomically Relevant Climate Indices via the Assemblage {VAR}},
  author  = {Barrette, Christophe and Goulet Coulombe, Philippe and Reinicke, Tim},
  journal = {The Econometrics Journal},
  year    = {2026},
  note    = {forthcoming}
}
```
