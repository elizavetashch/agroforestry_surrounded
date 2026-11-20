# Agroforestry yield analysis: targets pipeline
# Run with targets::tar_make(), inspect with targets::tar_visnetwork()

library(targets)

tar_option_set(packages = c("dplyr", "ggplot2", "mgcv"))
tar_source(here::here("02_scriptsandmodels", "functions"))

# settings -------------------------------------------------------------------

clim_vars <- c("temp_C_mean", "precip_mm_sum", "sun_MJ_m2_mean")
soil_vars <- c("clay", "sand", "silt")

zuur <- list(
  y      = "yield_tha",
  naive  = "yield_tha ~ distance_to_tree_strip + crop_unified",
  random = c("field", "year"),
  groups       = c("field", "year", "crop_unified", "crop_season", "distance_to_tree_strip"),
  covariates   = c("distance_to_tree_strip", "treeage", "AFage", "l_shdi", "l_contag", "l_ed",
                   "mean_slope", "sd_slope", "prop_swf_within", "swf_beyond", "PC1_c", "PC1_s"),
  condition    = "crop_season",
  interactions = c("crop_unified", "crop_season"),
  order_by     = c("field", "plot", "year")
)

infield_var <- "prop_swf_within"

# collinearity check (Zuur step 5) on the raw variables, before the PCA is built
collin_vars <- c("treeage", "AFage", "l_shdi", "l_contag", "l_ed", "mean_slope", "sd_slope",
                 "prop_swf_within", "swf_beyond", clim_vars, soil_vars)

family_formula <- "yield_tha ~ crop_unified + s(distance_to_tree_strip, k = 5, by = crop_unified) + s(field, bs = \"re\") + s(year, bs = \"re\")"

# the same exploration for the control fields, without the terms that need trees
zuur_ctrl <- modifyList(zuur, list(
  naive      = "yield_tha ~ crop_unified",
  groups     = setdiff(zuur$groups, "distance_to_tree_strip"),
  covariates = setdiff(zuur$covariates, c("distance_to_tree_strip", "treeage", "AFage"))
))
collin_vars_ctrl <- setdiff(collin_vars, c("treeage", "AFage"))
family_formula_ctrl <- "yield_tha ~ crop_unified + s(field, bs = \"re\") + s(year, bs = \"re\")"


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
  tar_target(f_fieldpolygons, path_in("fields", "fieldpolygons_with_control.csv"), format = "file"),
  tar_target(f_fertilization, path_in("fields", "DataDescription_Fertilization.csv"), format = "file"),
  tar_target(f_lulc, list.files(path_in("lulc"), "\\.tif$", recursive = TRUE, full.names = TRUE), format = "file"),
  tar_target(f_swf, list.files(path_in("swf_raw"), "\\.tif$", recursive = TRUE, full.names = TRUE), format = "file"),
  tar_target(f_soil, list.files(path_in("soilgrids"), "\\.tif$", full.names = TRUE), format = "file"),
  
  
  # (1) data preparation: single fields 
  tar_target(signal16, prep_signal16(f_signal16) |> check_nrow(120, "signal16")),
  tar_target(wendhausen, prep_wendhausen(signal16, f_wendhausen) |> check_field_data(140, "wendhausen")),
  tar_target(dornburg, prep_dornburg(signal16, f_dornburgvechta) |> check_field_data(280, "dornburg")),
  tar_target(vechta, prep_vechta(f_dornburgvechta) |> check_field_data(100, "vechta")),
  tar_target(reiffenhausen, prep_reiffenhausen(f_reiffenhausen) |> check_field_data(32, "reiffenhausen")),
  tar_target(mariensee, prep_mariensee(f_mariensee) |> check_field_data(60, "mariensee")),
  tar_target(forst, prep_forst(signal16, f_forst) |> check_field_data(80, "forst")),
  tar_target(gladbacherhof, prep_gladbacherhof(f_gladbacherhof) |> check_field_data(144, "gladbacherhof")),
  tar_target(koch, prep_koch(f_koch) |> check_field_data(451, "koch")),
  
  # (1) data preparation: merge and add site variables
  tar_target(merged,
             merge_fields(list(wendhausen, dornburg, vechta, reiffenhausen, mariensee, forst, gladbacherhof, koch),
                          f_fieldpolygons, f_fertilization)|>
               check_unique("id", "merged") |>
               check_no_na(c("yield_tha", "crop_unified", "latitude", "longitude", "treeage"), "merged") |>
               check_no_na(c("gs_start_date", "gs_end_date"), "merged", warn = TRUE)
  ),
  tar_target(merged_csv, save_table(merged, path_out("merged.csv")), format = "file"),
  tar_target(merged_with_control,
             merge_fields(list(wendhausen, dornburg, vechta, reiffenhausen, mariensee, forst, gladbacherhof, koch),
                          f_fieldpolygons, f_fertilization,
                          keep_control = TRUE) |>
               check_unique("id", "merged_with_control") |>
               check_no_na(c("yield_tha", "crop_unified", "latitude", "longitude", "treeage", "control_field"), "merged_with_control") |>
               check_no_na(c("gs_start_date", "gs_end_date"), "merged_with_control", warn = TRUE)
  ),
  tar_target(merged_with_control_csv, save_table(merged_with_control, path_out("merged_with_control.csv")), format = "file"),
  tar_target(field_coords, distinct(merged_with_control, field, control_field, latitude, longitude) |>
               check_unique(c("field", "control_field"), "field_coords")),
  tar_target(climate, extract_climate(field_coords, 2015:2023, path_in("dwd"))),
  tar_target(af_climate, add_climate(merged_with_control, climate) |>
               check_same_rows(merged_with_control, "af_climate") |>
               check_no_na(clim_vars, "af_climate", warn = TRUE)),
  tar_target(lulc, lulc_metrics(f_fieldpolygons, f_lulc) |> check_unique(c("field", "control_field", "year"), "lulc")),
  tar_target(af_lulc, add_lulc(af_climate, lulc) |>
               check_same_rows(af_climate, "af_lulc") |>
               check_no_na(c("l_shdi", "l_contag", "l_ed"), "af_lulc", warn = TRUE)),
  tar_target(af_slope, add_slope(af_lulc, f_fieldpolygons) |> check_same_rows(af_lulc, "af_slope")),
  tar_target(soil, soil_texture(f_fieldpolygons, f_soil) |> check_cols(soil_vars, "soil") |> check_unique(c("field", "control_field"), "soil")),
  tar_target(swf, swf_profiles(f_fieldpolygons, f_swf) |> check_unique(c("field", "control_field", "swf_year"), "swf")),
  tar_target(dffinal, finalise_dataset(af_slope, soil) |> check_unique("id", "dffinal")),
  tar_target(moddata_base, make_moddata(dffinal, swf, clim_vars, soil_vars, pca = FALSE)),
  tar_target(moddata_all, add_pca(moddata_base, clim_vars, soil_vars)),
  tar_target(moddata, select_system(moddata_all, "agroforestry") |> check_moddata()),
  tar_target(moddata_ctrl, select_system(moddata_all, "control") |> check_moddata("moddata_ctrl", control = TRUE)),
  tar_target(moddata_csv, save_table(moddata_all, path_out("moddata.csv")), format = "file"),
  tar_target(pca_files, pca_tables(moddata_all, clim_vars, soil_vars, res_dir("01_dataprep")), format = "file"),
  
  # (2) data exploration after Zuur et al. 2010 --------------------------------
  tar_target(zuur_dir, res_dir("02_exploration_zuur")),
  tar_target(zuur1, zuur_step1_outliers(moddata, zuur$y, zuur$covariates, zuur_dir), format = "file"),
  tar_target(zuur2, zuur_step2_variance(moddata, zuur$y, zuur$groups, zuur_dir), format = "file"),
  tar_target(zuur3, zuur_step3_normality(moddata, zuur$y, zuur$groups, zuur_dir), format = "file"),
  tar_target(zuur4, zuur_step4_zeros(moddata, zuur$y, zuur$groups, zuur_dir), format = "file"),
  tar_target(zuur5, zuur_step5_collinearity(moddata_base, collin_vars, zuur_dir), format = "file"),
  tar_target(zuur6, zuur_step6_relationships(moddata, zuur$y, zuur$covariates, zuur$condition, zuur_dir), format = "file"),
  tar_target(zuur7, zuur_step7_interactions(moddata, zuur$y, zuur$covariates, zuur$interactions, zuur_dir), format = "file"),
  tar_target(zuur8, zuur_step8_independence(moddata, zuur$naive, zuur$random, zuur$order_by, zuur_dir), format = "file"),
  tar_target(zuur9, family_qq(moddata, "yield_tha", zuur_dir), format = "file"),
  tar_target(family_check, family_compare(gamdata, family_formula, model_dir), format = "file"),
  tar_target(mean_sd_check, mean_sd_compare(moddata, zuur_dir), format = "file"),
  tar_target(add_season_croptype, season_croptype(moddata, zuur_dir), format = "file"), 
  tar_target(site_map, plot_site_map(moddata, dffinal, file.path(zuur_dir, "GermanMap.png")), format = "file"),
  tar_target(season_yield, plot_season_yield(moddata, zuur_dir), format = "file"),
  
  # (2b) data exploration after Zuur et al. 2010, control fields --------------------
  tar_target(zuur_ctrl_dir, res_dir("02_exploration_zuur_control")),
  tar_target(moddata_base_ctrl, select_system(moddata_base, "control")),
  tar_target(zuur1_ctrl, zuur_step1_outliers(moddata_ctrl, zuur_ctrl$y, zuur_ctrl$covariates, zuur_ctrl_dir), format = "file"),
  tar_target(zuur2_ctrl, zuur_step2_variance(moddata_ctrl, zuur_ctrl$y, zuur_ctrl$groups, zuur_ctrl_dir), format = "file"),
  tar_target(zuur3_ctrl, zuur_step3_normality(moddata_ctrl, zuur_ctrl$y, zuur_ctrl$groups, zuur_ctrl_dir), format = "file"),
  tar_target(zuur4_ctrl, zuur_step4_zeros(moddata_ctrl, zuur_ctrl$y, zuur_ctrl$groups, zuur_ctrl_dir), format = "file"),
  tar_target(zuur5_ctrl, zuur_step5_collinearity(moddata_base_ctrl, collin_vars_ctrl, zuur_ctrl_dir), format = "file"),
  tar_target(zuur6_ctrl, zuur_step6_relationships(moddata_ctrl, zuur_ctrl$y, zuur_ctrl$covariates, zuur_ctrl$condition, zuur_ctrl_dir), format = "file"),
  tar_target(zuur7_ctrl, zuur_step7_interactions(moddata_ctrl, zuur_ctrl$y, zuur_ctrl$covariates, zuur_ctrl$interactions, zuur_ctrl_dir), format = "file"),
  tar_target(zuur8_ctrl, zuur_step8_independence(moddata_ctrl, zuur_ctrl$naive, zuur_ctrl$random, zuur_ctrl$order_by, zuur_ctrl_dir), format = "file"),
  tar_target(zuur9_ctrl, family_qq(moddata_ctrl, "yield_tha", zuur_ctrl_dir), format = "file"),
  tar_target(family_check_ctrl, family_compare(gamdata_ctrl, family_formula_ctrl, zuur_ctrl_dir), format = "file"),
  tar_target(mean_sd_check_ctrl, mean_sd_compare(moddata_ctrl, zuur_ctrl_dir), format = "file"),
  tar_target(add_season_croptype_ctrl, season_croptype(moddata_ctrl, zuur_ctrl_dir), format = "file"),
  tar_target(season_yield_ctrl, plot_season_yield(moddata_ctrl, zuur_ctrl_dir), format = "file"),
  
  # (3) modelling ----------------------------------------------------------------
  tar_target(gamdata, add_swf_matrix(moddata)),
  tar_target(model_spec, model_formulas()),
  tar_target(models, fit_model(gamdata, model_spec), pattern = map(model_spec), iteration = "list"),
  
  # (4) model results visualisation -------------------------------------------
  tar_target(model_dir, res_dir("03_models")),
  tar_target(model_tables, save_model_tables(models, model_dir), format = "file"),
  tar_target(model_figures, save_model_figures(models, model_dir), format = "file"),
  tar_target(fig_relative_yield, plot_relative_yield(moddata, file.path(model_dir, "yield_distance.png")), format = "file"),
  tar_target(fig_aic, plot_aic(models, file.path(model_dir, "aic_h1_h2.png")), format = "file"),
  tar_target(fig_swf_function, plot_swf_function(models, gamdata, file.path(model_dir, "m01d_swf_function.png")), format = "file"),
  tar_target(fig_distance_crop, plot_distance_by_crop(models, gamdata, file.path(model_dir, "m01c_distance_by_crop.png")), format = "file"),
  tar_target(fig_h2_space_swf, plot_heterogeneity_space(moddata, "swf_beyond", infield_var,
                                                        file.path(model_dir, "h2_space_swf.png")), format = "file"),
  tar_target(fig_h2_space_shdi, plot_heterogeneity_space(moddata, "l_shdi", infield_var,
                                                         file.path(model_dir, "h2_space_shdi.png")), format = "file"),
  tar_target(fig_h2_raw_swf, plot_interaction_raw(moddata, "swf_beyond", infield_var,
                                                  file.path(model_dir, "h2_raw_swf.png")), format = "file"),
  tar_target(fig_h2_raw_shdi, plot_interaction_raw(moddata, "l_shdi", infield_var,
                                                   file.path(model_dir, "h2_raw_shdi.png")), format = "file"),
  tar_target(fig_h2_model_swf, plot_interaction_model(models, gamdata, "swf_beyond", infield_var,
                                                      file.path(model_dir, "h2_model_swf.png")), format = "file"),
  tar_target(fig_h2_model_shdi, plot_interaction_model(models, gamdata, "l_shdi", infield_var,
                                                       file.path(model_dir, "h2_model_shdi.png")), format = "file"),
  tar_target(fig_coef_forest, plot_coef_forest(models, file.path(model_dir, "effects_parametric_all_models.png")), format = "file"),
  tar_target(fig_smooth_heatmap, plot_smooth_heatmap(models, file.path(model_dir, "effects_smooth_all_models.png")), format = "file"),
  tar_target(fig_swf_profiles, plot_swf_profiles(swf, file.path(model_dir, "swf_profiles_by_field.png")), format = "file"),
  tar_target(fig_yield_by_field, plot_yield_by_field(moddata, file.path(model_dir, "yield_distance_by_field.png")), format = "file"),
  
  # (5) same models with relative yield (yield / maximum yield of the crop) ----
  tar_target(model_spec_rel, model_formulas("yield_rel_crop", crop_main = FALSE)),
  tar_target(models_rel, fit_model(gamdata, model_spec_rel), pattern = map(model_spec_rel), iteration = "list"),
  tar_target(model_dir_rel, res_dir("04_models_relative_yield")),
  tar_target(model_tables_rel, save_model_tables(models_rel, model_dir_rel), format = "file"),
  tar_target(model_figures_rel, save_model_figures(models_rel, model_dir_rel), format = "file"),
  tar_target(fig_aic_rel, plot_aic(models_rel, file.path(model_dir_rel, "aic_h1_h2.png")), format = "file"),
  tar_target(fig_coef_forest_rel, plot_coef_forest(models_rel, file.path(model_dir_rel, "effects_parametric_all_models.png")), format = "file"),
  tar_target(fig_smooth_heatmap_rel, plot_smooth_heatmap(models_rel, file.path(model_dir_rel, "effects_smooth_all_models.png")), format = "file"),
  tar_target(fig_relative_yield_rel, plot_relative_yield(moddata, file.path(model_dir_rel, "yield_distance.png"),
                                                         y = "yield_rel_crop", ylab = "Relative yield (yield / crop maximum)"), format = "file"),
  tar_target(fig_yield_by_field_rel, plot_yield_by_field(moddata, file.path(model_dir_rel, "yield_distance_by_field.png"),
                                                         y = "yield_rel_crop", ylab = "Relative yield (yield / crop maximum)"), format = "file"),
  tar_target(fig_distance_crop_rel, plot_distance_by_crop(models_rel, gamdata, file.path(model_dir_rel, "m01c_distance_by_crop.png"),
                                                          ylab = "Predicted relative yield"), format = "file"),
  tar_target(fig_swf_function_rel, plot_swf_function(models_rel, gamdata, file.path(model_dir_rel, "m01d_swf_function.png")), format = "file"),
  tar_target(fig_h2_model_swf_rel, plot_interaction_model(models_rel, gamdata, "swf_beyond", infield_var,
                                                          file.path(model_dir_rel, "h2_model_swf.png")), format = "file"),
  tar_target(fig_h2_model_shdi_rel, plot_interaction_model(models_rel, gamdata, "l_shdi", infield_var,
                                                           file.path(model_dir_rel, "h2_model_shdi.png")), format = "file"),
  
  
  # (6) control comparison 
  tar_target(control_dir, res_dir("05_control_comparison")),
  tar_target(c_yield, control_yield(merged_with_control)),
  
  tar_target(fig_yield_timeseries,
             plot_yield_timeseries(merged_with_control,
                                   file.path(control_dir, "yield_timeseries_af_control.png")),format = "file"),
  tar_target(fig_yield_bars,
             plot_yield_bars(merged_with_control,
                             file.path(control_dir, "yield_bars_af_control.png")),format = "file"),
  tar_target(fig_control_yield,
             plot_control_yield(c_yield,
                                 file.path(control_dir, "relative_yield_all_fields.png")),format = "file"),

  # (7) same models for the control fields --------------------------------------
  tar_target(gamdata_ctrl, add_swf_matrix(moddata_ctrl)),
  tar_target(model_spec_ctrl, model_formulas(control = TRUE)),
  tar_target(models_ctrl, fit_model(gamdata_ctrl, model_spec_ctrl), pattern = map(model_spec_ctrl), iteration = "list"),
  tar_target(model_dir_ctrl, res_dir("06_models_control")),
  tar_target(model_tables_ctrl, save_model_tables(models_ctrl, model_dir_ctrl), format = "file"),
  tar_target(model_figures_ctrl, save_model_figures(models_ctrl, model_dir_ctrl), format = "file"),
  tar_target(fig_aic_ctrl, plot_aic(models_ctrl, file.path(model_dir_ctrl, "aic_h1_h2.png")), format = "file"),
  tar_target(fig_coef_forest_ctrl, plot_coef_forest(models_ctrl, file.path(model_dir_ctrl, "effects_parametric_all_models.png")), format = "file"),
  tar_target(fig_smooth_heatmap_ctrl, plot_smooth_heatmap(models_ctrl, file.path(model_dir_ctrl, "effects_smooth_all_models.png")), format = "file"),
  tar_target(fig_swf_function_ctrl, plot_swf_function(models_ctrl, gamdata_ctrl, file.path(model_dir_ctrl, "m01d_swf_function.png")), format = "file"),
  tar_target(fig_h2_model_swf_ctrl, plot_interaction_model(models_ctrl, gamdata_ctrl, "swf_beyond", infield_var,
                                                           file.path(model_dir_ctrl, "h2_model_swf.png")), format = "file"),
  tar_target(fig_h2_model_shdi_ctrl, plot_interaction_model(models_ctrl, gamdata_ctrl, "l_shdi", infield_var,
                                                            file.path(model_dir_ctrl, "h2_model_shdi.png")), format = "file"),
  
  # (7b) relative yield models for the control fields ------------------------------
  tar_target(model_spec_ctrl_rel, model_formulas("yield_rel_crop", crop_main = FALSE, control = TRUE)),
  tar_target(models_ctrl_rel, fit_model(gamdata_ctrl, model_spec_ctrl_rel), pattern = map(model_spec_ctrl_rel), iteration = "list"),
  tar_target(model_dir_ctrl_rel, res_dir("07_models_control_relative_yield")),
  tar_target(model_tables_ctrl_rel, save_model_tables(models_ctrl_rel, model_dir_ctrl_rel), format = "file"),
  tar_target(model_figures_ctrl_rel, save_model_figures(models_ctrl_rel, model_dir_ctrl_rel), format = "file"),
  tar_target(fig_coef_forest_ctrl_rel, plot_coef_forest(models_ctrl_rel, file.path(model_dir_ctrl_rel, "effects_parametric_all_models.png")), format = "file"),
  tar_target(fig_smooth_heatmap_ctrl_rel, plot_smooth_heatmap(models_ctrl_rel, file.path(model_dir_ctrl_rel, "effects_smooth_all_models.png")), format = "file"),

  # (8) agroforestry vs control: field profiles and model effects ---------------
  tar_target(fig_field_profiles, plot_field_profiles(moddata_all, file.path(control_dir, "field_profiles_yield.png")), format = "file"),
  tar_target(fig_field_profiles_rel, plot_field_profiles(moddata_all, file.path(control_dir, "field_profiles_relative_yield.png"),
                                                         y = "yield_rel_crop", ylab = "Relative yield (yield / crop maximum)"), format = "file"),
  tar_target(fig_effects_compare, plot_effects_compare(models, models_ctrl, file.path(control_dir, "effects_af_vs_control_m02d.png")), format = "file"),
  tar_target(fig_coef_forest_sets, plot_coef_forest_sets(list("Yield (t/ha)" = models, "Relative yield" = models_rel,
                                                              "Control" = models_ctrl, "Control relative yield" = models_ctrl_rel),
                                                         file.path(control_dir, "effects_parametric_all_sets.png"),
                                                         colours = c("#00A087", "#3C5488", "#7E6148", "#E64B35")), format = "file")
  )