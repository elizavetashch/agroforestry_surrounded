# Agroforestry yield: in-field and beyond

Targets approach based on the lectures given by Ondřej Mottl.

## Folder structure

```
_targets.R                    pipeline definition and settings
01_inputdata/                 raw input data (never changed by the code)
02_scriptsandmodels/
  functions/                  all functions used by the pipeline
03_outputdata/                moddata.csv (final model dataset)
04_resultsfigurestables/
  01_dataprep/                PCA loadings and explained variance
  02_exploration_zuur/        Zuur steps 1-8: figures and tables
  03_models/                  AIC table, model terms, k-check, diagnostics, smooth plots
```

## Run

```r
targets::tar_manifest()     # list of all steps
targets::tar_visnetwork()   # dependency graph
targets::tar_make()         # run everything that is outdated
targets::tar_read(moddata)  # load a result
```


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
