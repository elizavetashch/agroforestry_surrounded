# Agroforestry yield analysis: targets pipeline
# Run with targets::tar_make(), inspect with targets::tar_visnetwork()

library(targets)

tar_option_set(packages = c("dplyr", "ggplot2", "mgcv"))
tar_source(here::here("02_scriptsandmodels", "functions"))

# settings -------------------------------------------------------------------

clim_vars <- c("temp_C_mean", "precip_mm_sum", "sun_MJ_m2_mean")
soil_vars <- c("clay", "sand", "silt")

zuur <- list(
  y            = "yield_rel",
  groups       = c("field", "year", "crop_unified", "crop_season", "distance_to_tree_strip"),
  covariates   = c("distance_to_tree_strip", "treeage", "AFage", "l_shdi", "l_contag", "l_ed",
                   "mean_slope", "prop_swf_within", "PC1_c", "PC1_s"),
  condition    = "crop_season",
  interactions = c("crop_unified", "crop_season"),
  naive        = "yield_rel ~ distance_to_tree_strip",
  random       = c("field", "year", "crop_season"),
  order_by     = c("field", "plot", "year")
)

list(
  # (1) data preparation: input files ------------------------------------------
  tar_target(f_signal16, path_in("yields", "BONARES_SIGNAL16",
                                 "signal.ID_7048_BIOMASSES_SIGNAL_PROJECT_V1_APR_08_2020.csv"), format = "file"),
  tar_target(f_wendhausen, c(
    path_in("yields", "BONARES_Cropland agroforestry 2019-2020", "signal.ID_7069_DATEN_WENDHAUSEN_2019_20.csv"),
    path_in("yields", "BONARES_Cropland agroforestry 2021", "signal.ID_7070_DATEN_WENDHAUSEN_2021.csv"),
    path_in("yields", "BONARES_Cropland agroforestry 2022", "signal.ID_7071_DATEN_WENDHAUSEN_2022.csv"),
    path_in("yields", "BONARES_Cropland agroforestry 2023", "signal.ID_7077_DATEN_WENDHAUSEN_2023_CR.csv")
  ), format = "file"),
  tar_target(f_dornburgvechta, path_in("yields", "BONARES_DornburgVechta1823", "signal.ID_7088_CROP_YIELD.csv"), format = "file"),
  tar_target(f_reiffenhausen, path_in("yields", "BONARES_Reiffenhausen16", "signal.ID_7039_REIFFENHAUSEN_BIOMASS_DATA_V2.csv"), format = "file"),
  tar_target(f_mariensee, path_in("yields", "Mariensee1719", "signal.ID_7041_BIOMASSE_17_18_MS_280319.csv"), format = "file"),
  tar_target(f_forst, path_in("yields", "BONARES_Forst1920", "signal.ID_7060_CROP_YIELDS_FORST_2019_2020.csv"), format = "file"),
  tar_target(f_gladbacherhof, path_in("yields", "Gladbacherhof21", "yield.csv"), format = "file"),
  tar_target(f_koch, path_in("yields", "Koch25", "yields_wintercrops.csv"), format = "file"),
  tar_target(f_fieldpolygons, path_in("fields", "fieldpolygons.csv"), format = "file"),
  tar_target(f_fertilization, path_in("fields", "DataDescription_Fertilization.csv"), format = "file"),
  tar_target(f_slope, path_in("slope", "slope.csv"), format = "file"),
  tar_target(f_lulc, list.files(path_in("lulc"), "\\.tif$", recursive = TRUE, full.names = TRUE), format = "file"),
  tar_target(f_swf, list.files(path_in("swf"), "\\.tif$", recursive = TRUE, full.names = TRUE), format = "file"),
  tar_target(f_soil, list.files(path_in("soilgrids"), "\\.tif$", full.names = TRUE), format = "file"),

  # (1) data preparation: single fields (row numbers from the original scripts)
  tar_target(signal16, prep_signal16(f_signal16) |> check_nrow(120, "signal16")),
  tar_target(wendhausen, prep_wendhausen(signal16, f_wendhausen) |> check_field_data(140, "wendhausen")),
  tar_target(dornburg, prep_dornburg(signal16, f_dornburgvechta) |> check_field_data(280, "dornburg")),
  tar_target(vechta, prep_vechta(f_dornburgvechta) |> check_field_data(100, "vechta")),
  tar_target(reiffenhausen, prep_reiffenhausen(f_reiffenhausen) |> check_field_data(32, "reiffenhausen")),
  tar_target(mariensee, prep_mariensee(f_mariensee) |> check_field_data(60, "mariensee")),
  tar_target(forst, prep_forst(signal16, f_forst) |> check_field_data(80, "forst")),
  tar_target(gladbacherhof, prep_gladbacherhof(f_gladbacherhof) |> check_field_data(144, "gladbacherhof")),
  tar_target(koch, prep_koch(f_koch) |> check_field_data(480, "koch")),

  # (1) data preparation: merge and add site variables
  tar_target(merged,
    merge_fields(list(wendhausen, dornburg, vechta, reiffenhausen, mariensee, forst, gladbacherhof, koch),
                 f_fieldpolygons, f_fertilization) |>
      check_unique("id", "merged") |>
      check_no_na(c("yield_tha", "crop_unified", "latitude", "longitude", "treeage"), "merged") |>
      check_no_na(c("gs_start_date", "gs_end_date"), "merged", warn = TRUE)
  ),
  tar_target(field_coords, distinct(merged, field, latitude, longitude) |> check_unique("field", "field_coords")),
  tar_target(climate, extract_climate(field_coords, 2015:2023, path_in("dwd"))),
  tar_target(af_climate, add_climate(merged, climate) |>
               check_same_rows(merged, "af_climate") |>
               check_no_na(clim_vars, "af_climate", warn = TRUE)),
  tar_target(lulc, lulc_metrics(f_fieldpolygons, f_lulc) |> check_unique(c("field", "year"), "lulc")),
  tar_target(af_lulc, add_lulc(af_climate, lulc) |>
               check_same_rows(af_climate, "af_lulc") |>
               check_no_na(c("l_shdi", "l_contag", "l_ed"), "af_lulc", warn = TRUE)),
  tar_target(af_slope, add_slope(af_lulc, f_slope) |> check_same_rows(af_lulc, "af_slope")),
  tar_target(soil, soil_texture(f_fieldpolygons, f_soil) |> check_cols(soil_vars, "soil") |> check_unique("field", "soil")),
  tar_target(swf, swf_profiles(f_fieldpolygons, f_swf) |> check_unique(c("field", "swf_year"), "swf")),
  tar_target(dffinal, finalise_dataset(af_slope, soil) |> check_unique("id", "dffinal")),
  tar_target(moddata, make_moddata(dffinal, swf, clim_vars, soil_vars) |> check_moddata()),
  tar_target(moddata_csv, save_table(moddata, path_out("moddata.csv")), format = "file"),
  tar_target(pca_files, pca_tables(moddata, clim_vars, soil_vars, res_dir("01_dataprep")), format = "file"),

  # (2) data exploration after Zuur et al. 2010 --------------------------------
  tar_target(zuur_dir, res_dir("02_exploration_zuur")),
  tar_target(zuur1, zuur_step1_outliers(moddata, zuur$y, zuur$covariates, zuur_dir), format = "file"),
  tar_target(zuur2, zuur_step2_variance(moddata, zuur$y, zuur$groups, zuur_dir), format = "file"),
  tar_target(zuur3, zuur_step3_normality(moddata, zuur$y, zuur$groups, zuur_dir), format = "file"),
  tar_target(zuur4, zuur_step4_zeros(moddata, zuur$y, zuur$groups, zuur_dir), format = "file"),
  tar_target(zuur5, zuur_step5_collinearity(moddata, zuur$covariates, zuur_dir), format = "file"),
  tar_target(zuur6, zuur_step6_relationships(moddata, zuur$y, zuur$covariates, zuur$condition, zuur_dir), format = "file"),
  tar_target(zuur7, zuur_step7_interactions(moddata, zuur$y, zuur$covariates, zuur$interactions, zuur_dir), format = "file"),
  tar_target(zuur8, zuur_step8_independence(moddata, zuur$naive, zuur$random, zuur$order_by, zuur_dir), format = "file"),

  # (3) modelling ----------------------------------------------------------------
  tar_target(gamdata, add_swf_matrix(moddata)),
  tar_target(model_spec, model_formulas()),
  tar_target(models, fit_model(gamdata, model_spec), pattern = map(model_spec), iteration = "list"),

  # (4) model results visualisation -------------------------------------------
  tar_target(model_dir, res_dir("03_models")),
  tar_target(model_tables, save_model_tables(models, model_dir), format = "file"),
  tar_target(model_figures, save_model_figures(models, model_dir), format = "file"),
  tar_target(fig_relative_yield, plot_relative_yield(moddata, file.path(model_dir, "relative_yield_distance.png")), format = "file"),
  tar_target(fig_distance_crop, plot_distance_by_crop(models, gamdata, file.path(model_dir, "m_01_distance_by_crop.png")), format = "file")
)
