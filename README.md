# Agroforestry yield: in-field and beyond

Master thesis project (Elizaveta Shcherbinina). The whole workflow is a
[`targets`](https://docs.ropensci.org/targets/) pipeline: data preparation,
data exploration after Zuur et al. (2010), modelling, and model visualisation.

## Folder structure

```
_targets.R                    pipeline definition and settings
01_inputdata/                 raw input data (never changed by the code)
02_scriptsandmodels/
  functions/                  all functions used by the pipeline
  tests/                      unit tests and validation against the old pipeline
03_outputdata/                moddata.csv (final model dataset)
04_resultsfigurestables/
  01_dataprep/                PCA loadings and explained variance
  02_exploration_zuur/        Zuur steps 1-8: figures and tables
  03_models/                  AIC table, model terms, k-check, diagnostics, smooth plots
```

Intermediate datasets are not written to disk. They are kept in the `_targets/`
store and can be loaded with `targets::tar_read(name)`, for example
`tar_read(merged)` or `tar_read(dffinal)`.

## Setup

1. Open `MA_agroforestry_targets.Rproj` in RStudio (this sets the project root for `here::here()`).
2. Install the packages once:
   ```r
   install.packages(c("targets", "here", "dplyr", "tidyr", "readr", "stringr", "janitor",
                      "ggplot2", "mgcv", "terra", "sf", "landscapemetrics", "rdwd", "visNetwork"))
   ```
3. Copy the raw data into `01_inputdata/` as listed in `01_inputdata/README.md`.

## Run

```r
targets::tar_manifest()     # list of all steps
targets::tar_visnetwork()   # dependency graph
targets::tar_make()         # run everything that is outdated
targets::tar_read(moddata)  # load a result
```

Only steps whose code or input changed are rerun. The first run downloads the
DWD climate grids into `01_inputdata/dwd/` (about 10 min); later runs reuse them.

## Pipeline

| Part | Targets | Function file |
|---|---|---|
| (1) Data preparation | `signal16`, `wendhausen`, ... `koch` | `01_prep_fields.R` |
| | `merged` | `02_prep_merge.R` |
| | `climate`, `af_climate`, `lulc`, `af_lulc`, `af_slope`, `soil`, `swf` | `03_prep_spatial.R` |
| | `dffinal`, `moddata` → `03_outputdata/moddata.csv` | `04_prep_moddata.R` |
| (2) Exploration (Zuur) | `zuur1` ... `zuur8` | `05_exploration_zuur.R` |
| (3) Modelling | `gamdata`, `model_spec`, `models` | `06_models.R` |
| (4) Visualisation | `model_tables`, `model_figures`, `fig_relative_yield`, `fig_distance_crop` | `07_model_visualisation.R` |

Settings (climate and soil variables, Zuur variables) are at the top of `_targets.R`.
Model formulas are in `model_formulas()` in `06_models.R`.

## Tests

* **Checks inside the pipeline** (`functions/00_checks.R`): every new dataset is
  checked when it is created (expected row numbers per field, required columns,
  unique ids, no lost or duplicated rows after joins, value ranges, missing
  values). A failing check stops `tar_make()` and names the dataset and the problem.
  Warnings (e.g. missing covariates) are listed with
  `targets::tar_meta(fields = warnings, complete_only = TRUE)`.
* **Unit tests** of the helper functions on small made-up data:
  `source(here::here("02_scriptsandmodels", "tests", "test_functions.R"))`
* **Validation against the old pipeline**: after `tar_make()`, set `old_project`
  at the top of `02_scriptsandmodels/tests/validate_against_old.R` and source it.
  It compares `merged`, `af_climate`, `af_lulc`, `af_slope`, `dffinal` and `moddata`
  with `AF_merged.csv`, `AF_climate.csv`, `AF_lulc.csv`, `AF_slope.csv`,
  `dffinal_*.csv` and `moddata.csv` (rows, columns, values) and saves a table
  to `04_resultsfigurestables/00_validation/`.

## Reference

Zuur AF, Ieno EN, Elphick CS (2010) A protocol for data exploration to avoid
common statistical problems. Methods in Ecology and Evolution 1: 3-14.
https://doi.org/10.1111/j.2041-210X.2009.00001.x
