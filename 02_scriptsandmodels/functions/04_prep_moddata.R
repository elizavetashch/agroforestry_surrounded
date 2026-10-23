af_planting_year <- c(
  Dornburg = 2007, Forst = 2010, Gladbacherhof = 2020, IhingerHof = 2008,
  Mariensee = 2008, Reiffenhausen = 2011, Vechta = 2019, Wendhausen = 2008
)

finalise_dataset <- function(df, soil) {
  d <- dplyr::left_join(df, soil, by = "field", relationship = "many-to-one")
  fert_raw <- d$fertilization_rate_kg_n_p_k_ha_1_year_1
  for (dash in c("\u2013", "\u2212")) fert_raw <- gsub(dash, "-", fert_raw, fixed = TRUE)
  d$fert_N <- suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", fert_raw)))
  d$fert_N[fert_raw == "no fertilization"] <- 0
  d$year_AFplanting <- unname(af_planting_year[d$field])
  d$AFage <- d$year - d$year_AFplanting
  d <- unique(d)
  d$ID_nodist <- gsub("_\\d+(?:\\.\\d+)?m_", "_", d$id, perl = TRUE)
  d
}

swf_reference_year <- function(year) {
  dplyr::case_when(
    year == 2016         ~ 2015,
    year %in% 2017:2019  ~ 2018,
    year %in% 2020:2024  ~ 2021
  )
}

first_pc <- function(units, vars) {
  stats::prcomp(units[vars], scale. = TRUE)
}

make_moddata <- function(df, swf, clim_vars, soil_vars) {
  keys140 <- c("field", "year", "crop_unified", "distance_to_tree_strip")
  const_vars <- c("AFage", "treeage", "l_shdi", "l_ed", "l_np", "l_contag", "mean_slope",
                  "prop_swf_within", "sowing_month", "harvest_month")
  swf_cols <- grep("^swf_d", names(swf), value = TRUE)

  d <- df |>
    dplyr::mutate(swf_year = swf_reference_year(year)) |>
    dplyr::left_join(swf, by = c("field", "swf_year"), relationship = "many-to-one") |>
    dplyr::select(plot, dplyr::all_of(c(keys140, clim_vars, soil_vars, const_vars)),
                  swf_year, dplyr::all_of(swf_cols), yield_tha) |>
    dplyr::distinct() |>
    dplyr::mutate(
      crop_season = dplyr::case_when(
        sowing_month %in% c("Sep", "Oct", "Nov", "Dec", "Jan") ~ "winter",
        sowing_month %in% c("Mar", "Apr", "May")               ~ "summer"
      ),
      key34      = paste(field, year, crop_unified, sep = "_"),
      field_year = paste(field, year, sep = "_"),
      key140     = paste(key34, distance_to_tree_strip, sep = "_")
    ) |>
    dplyr::mutate(key140_mean = mean(yield_tha), n_cell = dplyr::n(), .by = key140)

  key34_ref <- d |>
    dplyr::distinct(key34, key140, key140_mean) |>
    dplyr::summarise(key34_mean = mean(key140_mean), .by = key34)

  d <- d |>
    dplyr::left_join(key34_ref, by = "key34") |>
    dplyr::mutate(yield_rel = yield_tha / key34_mean, key140_rel = key140_mean / key34_mean) |>
    dplyr::mutate(dist_rel = distance_to_tree_strip / max(distance_to_tree_strip), .by = field_year)

  clim_units <- dplyr::distinct(d, dplyr::across(dplyr::all_of(c("key34", clim_vars))))
  soil_units <- dplyr::distinct(d, dplyr::across(dplyr::all_of(c("key34", soil_vars))))
  clim_units$PC1_c <- first_pc(clim_units, clim_vars)$x[, 1]
  soil_units$PC1_s <- first_pc(soil_units, soil_vars)$x[, 1]

  n <- nrow(d)
  d <- d |>
    dplyr::left_join(dplyr::select(clim_units, key34, PC1_c), by = "key34") |>
    dplyr::left_join(dplyr::select(soil_units, key34, PC1_s), by = "key34") |>
    dplyr::mutate(dplyr::across(c(field, year, crop_unified, crop_season), as.factor))
  if (nrow(d) != n) stop("make_moddata: PCA join duplicated rows (more than one climate or soil value per field x year x crop)")
  d
}

pca_tables <- function(moddata, clim_vars, soil_vars, dir) {
  one <- function(vars, label) {
    units <- dplyr::distinct(moddata, dplyr::across(dplyr::all_of(c("key34", vars))))
    p <- first_pc(units, vars)
    imp <- summary(p)$importance
    load <- data.frame(pca = label, variable = rownames(p$rotation), p$rotation, row.names = NULL)
    var <- data.frame(pca = label, statistic = rownames(imp), imp, row.names = NULL)
    list(load, var)
  }
  c1 <- one(clim_vars, "climate")
  s1 <- one(soil_vars, "soil")
  c(save_table(dplyr::bind_rows(c1[[1]], s1[[1]]), file.path(dir, "pca_loadings.csv")),
    save_table(dplyr::bind_rows(c1[[2]], s1[[2]]), file.path(dir, "pca_variance.csv")))
}

# matrix columns needed for the functional SWF term s(R, by = SWF)
add_swf_matrix <- function(moddata) {
  d <- as.data.frame(moddata)
  swf_cols <- grep("^swf_d\\d+$", names(d), value = TRUE)
  radii <- as.numeric(sub("swf_d", "", swf_cols))
  swf_cols <- swf_cols[order(radii)]
  radii <- sort(radii)
  d$SWF <- as.matrix(d[swf_cols])
  d$R <- matrix(radii, nrow = nrow(d), ncol = length(radii), byrow = TRUE)
  d[swf_cols] <- NULL
  d
}
