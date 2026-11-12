# Data dictionary

All files in this folder (not in `data/raw/`) are comma-separated text with a header row. Dates are the first day of the month
(`YYYY-MM-01`). Empty cells are missing values. The panels are the cleaned inputs used by the
estimator. Their provenance and transformations are documented below so that the supplied analysis
inputs can be audited. The package loader (`code/funcs/load_data.R`) re-derives the transformed panels from the
level panels at run time and stops if they differ from the shipped `*_transformed.csv` files.

The optional step `code/00_prepare_data.R` rebuilds all files from the raw source files in `data/raw/`
(see [Data preparation](#data-preparation) below). It is not part of the main run, and its packages
are not needed to run this package.

Two raw source files, the ACI workbook and the Macrobond export of the NOAA series, were obtained under
licences that do not permit redistribution. `data/raw/README.md` lists every raw file with its source
and explains how to obtain each one from its provider. The seasonally adjusted and transformed data
used in the paper are included here for replication purposes and are sufficient to reproduce all reported results.

## Sources

| Source | Series | Vintage / access | Licence |
|---|---|---|---|
| Actuaries Climate Index (ACI), U.S. portion, monthly, unstandardised components | national: Sea.Level, CDD, Rx5Day, T10, T90, WP90. Regional: Sea.Level, Rx5Day, T10, T90, WP90 for the regions ALA, CEA, CWP, MID, SEA, SPL, SWP (no sea level for MID) | downloaded from https://actuariesclimateindex.org (data downloads, "unstandardized" components), January 2026; workbook version 1.1, values calculated 2025-11-17 | ACI terms of use (https://actuariesclimateindex.org/terms-of-use/) |
| NOAA National Centers for Environmental Information, Climate at a Glance (regional time series), and NOAA Global Monitoring Laboratory (Mauna Loa CO2) | Temperature, Temperature average, Precipitation, Heating and Cooling Degree Days, Palmer indices (PDSI, PHDI, PMDI, Z-Index) for the nine U.S. climate regions and the contiguous United States. CO2 average and trend | monthly, 1961-01 to 2025-12, retrieved via a Macrobond export, February 2026 | U.S. Government work (public domain). The Macrobond export itself is licensed and not redistributed |
| FRED-MD (McCracken and Ng), monthly database | INDPRO, PCEPI, UNRATE, EXUSUKx, HOUST, and FEDFUNDS (for SHADOWRATE outside 1990-04 to 2022-02) | vintage 2025-12 (`data/raw/2025-12-md.csv`, as published) | Federal Reserve Bank of St. Louis, public |
| Wu-Xia shadow federal funds rate | SHADOWRATE, 1990-04 to 2022-02 | Federal Reserve Bank of Atlanta, monthly (`data/raw/WuXiaShadowRate.xlsx`, as published) | public |
| NBER business cycle dates | recession start and end months | https://www.nber.org/research/data/us-business-cycle-expansions-and-contractions, accessed 1 October 2026 | public |

## Files

### `aci_national_level.csv` (773 months, 1961-01 to 2025-05)
Columns `Sea.Level`, `CDD`, `Rx5Day`, `T10`, `T90`, `WP90`: the six ACI national components,
seasonally adjusted with X-13ARIMA-SEATS (`RJDemetra::x13`, default specification, series treated as
monthly from 1961-01). CDD is consecutive dry days (reported annually and interpolated to monthly
frequency by the ACI). It is used only for the headline ACI benchmark and is dropped from the
constituent panel (see paper, footnote in Section 3.1).

### `aci_national_transformed.csv` (772 months, 1961-02 to 2025-05)
Columns `Sea.Level`, `Rx5Day`, `T10`, `T90`, `WP90`. Transformation, for each column x:
1. month-on-month log difference in percent, `100 * (log x_t - log x_{t-1})`.
2. standardisation with the monthly mean and standard deviation of step 1 over the reference years
   1961-1990 (the ACI normalisation convention): `(y_t - mean_m) / sd_m` where m is the calendar month.
3. `T10` is multiplied by -1 so that all components increase with adverse climate stress.
4. `CDD` is dropped.

### `aci_headline.csv` (772 months)
Column `ACI`: the equal-weight benchmark used in the VAR for all specifications, i.e. the row mean of
the five columns of `aci_national_transformed.csv` (not five-year smoothed, unlike the published ACI).

### `aci_regional_level.csv` (773 months) and `aci_regional_transformed.csv` (772 months)
34 columns named `<component>.<region>` (regions: ALA Alaska, CEA Central East Atlantic, CWP Central
West Pacific, MID Midwest, SEA Southeast Atlantic, SPL Southern Plains, SWP Southwest Pacific).
Level file: seasonally adjusted as above. Transformed file: absolute value of the level, then steps
1-3 above (log difference, 1961-1990 monthly standardisation, sign flip of all `T10.*` columns). The
absolute value is needed because some seasonally adjusted regional anomalies are negative.

### `noaa_level.csv` (780 months, 1961-01 to 2025-12) and `noaa_transformed.csv` (779 months)
82 columns with the original NOAA/Macrobond series names. Level file: seasonally adjusted with
X-13ARIMA-SEATS where the series allowed it (`RJDemetra::x13`). Series with missing values, fewer than
five years of data, or zero variance are left unadjusted. All 82 series in the panel were adjusted. Internal gaps of at most two months were
linearly interpolated before adjustment. Transformed file, for each column:
1. month-on-month change: log difference for series that are strictly positive, first difference
   otherwise. The rule per series is in `noaa_transformation_rules.csv` (35 log differences, 47 first
   differences). No scaling by 100.
2. the first month is dropped. Non-finite values are set to missing and internal gaps linearly
   interpolated. Series with more than 7.5% missing values or zero variance were excluded (this
   selection produced the 82 series).
3. monthly standardisation on 1961-1990 as for the ACI. Non-finite standardised values are set to 0.

### `noaa_transformation_rules.csv`
`name`, `rule` in {`logdiff`, `diff`}: the transformation applied to each NOAA series.

### `macro_fredmd_2025-12.csv` (802 months, 1959-01 to 2025-10) and `macro_transformations.csv`
Levels of `INDPRO`, `PCEPI`, `UNRATE`, `EXUSUKx`, `HOUST` as in FRED-MD vintage 2025-12 (its last month,
2025-11, is incomplete and dropped). `SHADOWRATE` is the Wu-Xia shadow federal funds rate for 1990-04 to
2022-02 and the FRED-MD effective federal funds rate `FEDFUNDS` before and after. Transformations
applied by the loader: `100 * dlog` for INDPRO, PCEPI, EXUSUKx, HOUST. First difference for UNRATE.
SHADOWRATE in levels. FRED-MD transformation codes are listed for reference in
`macro_transformations.csv` (for SHADOWRATE the code of FEDFUNDS).

### Label files
- `labels_aci_national.csv`: `index`, `name` (position of each national component in the weight vector).
- `labels_aci_regional.csv`: `index`, `name`, `variable`, `region`.
- `labels_noaa.csv`: `index`, `name`, `region_raw`, `variable_raw`, `region_agg`, `variable_agg`. The
  aggregated columns implement Tables C.1 and C.2 of the paper: regions Upper Midwest (Up.MidW.), Ohio
  Valley (Ohio.V.), Northwest (N.West), Northern Rockies & Plains (N.Rockies), Northeast (N.East) are
  kept. All other geographic categories are `Rest`. Variable families Temperature (Temp., including
  "Temperature, Average"), Precipitation (Precip.), Heating Degree Days (Heat.DD), Palmer Z-Index
  (Palmer.Z) are kept. CO2 level and trend, Cooling Degree Days and the three other Palmer indices are
  `Others`.

### `nber_recessions.csv`
`start`, `end`: NBER peak and trough month of each of the nine recessions between 1960 and 2020.

## Estimation sample
1975-01 to 2025-05 on the calendar. The months 2020-01 to 2021-06 (COVID-19) are removed. The last
month (2025-05) is held out, so estimation uses observations up to 2025-04 (587 months in the sample,
586 used, 12 lags, 561 gap-valid observation pairs for the likelihood criterion).

## Data preparation
`code/00_prepare_data.R` documents how the files above were built from the raw source files. It is
optional and not called by `run_all.R`. Run it from the package root:

    Rscript code/00_prepare_data.R

It reads the raw files in `data/raw/` (listed, with sources and licences, in `data/raw/README.md`, and
with sizes and SHA-256 checksums in `data/raw/MANIFEST.csv`; the ACI workbook and the Macrobond export
must first be obtained and placed there; the FRED-MD vintage and the Wu-Xia workbook are included as
published), seasonally adjusts the ACI and NOAA series with X-13ARIMA-SEATS, selects and transforms the
series as described above (the transformed panels with the loader's own functions), builds `SHADOWRATE`
from the Wu-Xia shadow rate and FRED-MD `FEDFUNDS`, and writes all files to `data/_prepared/`. The
shipped files in `data/` are never overwritten. It then compares every rebuilt file with the shipped
one and prints the maximal absolute and relative differences. The files that do not depend on the
seasonal adjustment (series selection and transformation rules, labels, macro data, recession dates)
must match, otherwise the script stops. X-13ARIMA-SEATS (JDemetra+, run in Java) is deterministic on a
given machine, but its results can differ slightly across platforms and Java versions. Differences in
the seasonally adjusted and transformed panels are therefore reported per file, with the number of
series affected, rather than treated as errors. It needs the optional packages `RJDemetra`, `rJava` and
`readxl` (see `requirements.csv`) and a Java runtime, version 8 or later, for `RJDemetra`. It takes
5-10 minutes on a laptop (AMD Ryzen 7 PRO 5850U, 32 GB RAM).
