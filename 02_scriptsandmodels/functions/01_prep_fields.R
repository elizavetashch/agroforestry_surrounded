# Cleaning of the raw yield datasets, one function per data source.
# Column positions (e.g. d[, c(1:5, 9:12)]) follow the raw files exactly.

read_raw <- function(file, data_id = NULL, yield_unit = NULL) {
  d <- readr::read_csv(file, show_col_types = FALSE)
  if (!is.null(data_id)) d$data_id <- data_id
  if (!is.null(yield_unit)) d$yield_unit <- yield_unit
  janitor::clean_names(d)
}

drop_na_yield <- function(d) d[!is.na(d$yield), ]

# SIGNAL 2016-2017: Wendhausen, Dornburg, Forst (Mariensee is taken from the 2017-19 file)
prep_signal16 <- function(file) {
  d <- read_raw(file, data_id = "signal16")
  d <- rename_col(d, "harvest_year", "year")
  d$yield_unit <- "g/m2 per year"
  d <- rename_col(d, "product_crop_grass_row_and_mono", "crop")
  d <- rename_col(d, "straw_or_grass_dry_mass_2016", "yield_straw")
  d <- rename_col(d, "grain_corn_or_grass_dry_mass_16", "yield")
  d <- d[, c(1:5, 9:12)]
  d$distance_to_tree_strip <- dist_from_id(d$id)
  d <- d[!(is.na(d$yield_straw) & is.na(d$yield)), ]
  d$field <- dplyr::case_when(
    d$site == "Dornburg"   ~ "Dornburg",
    d$site == "Forst"      ~ "Forst",
    d$site == "Wendhausen" ~ "Wendhausen"
  )
  d <- d[!is.na(d$field), ]
  rename_col(d, "id", "plot")
}

prep_wendhausen <- function(signal16, files) {
  w1920 <- read_raw(pick_file(files, "2019_20"), "wendhausen1920", "t/ha")
  w1920 <- w1920[, -(10:15)]
  w1920$crop <- "silage maize"
  w1920 <- drop_na_yield(rename_col(w1920, "sm_dm", "yield"))

  w21 <- read_raw(pick_file(files, "_2021"), "wendhausen21", "t/ha")
  w21 <- rename_col(w21, "sb_dm_straw", "yield_straw")
  w21 <- rename_col(w21, "sb_dm", "yield")
  w21$crop <- "summer barley"
  w21 <- drop_na_yield(w21[, -(11:15)])

  w22 <- read_raw(pick_file(files, "_2022"), "wendhausen22", "t/ha")
  w22 <- w22[, -(10:14)]
  w22$crop <- "oil rape"
  w22 <- drop_na_yield(rename_col(w22, "or_dm", "yield"))

  w23 <- read_raw(pick_file(files, "_2023"), "wendhausen23", "t/ha")
  w23 <- w23[, -(10:14)]
  w23$crop <- "sillage maize"
  w23 <- drop_na_yield(rename_col(w23, "sm_dm", "yield"))

  sig <- signal16[signal16$field == "Wendhausen", ]
  sig$lat <- 52.33336
  sig$lon <- 10.63236

  d <- dplyr::bind_rows(sig, w1920, w21, w22, w23)
  d$tree_species <- "poplar"
  d$field <- "Wendhausen"
  d
}

prep_dornburg <- function(signal16, file) {
  d <- read_raw(file)
  d$distance_to_tree_strip <- dist_from_id(d$id)
  d$data_id <- "dornburg1823"
  d$tree_species <- "poplar"
  d$yield_unit <- "g/m2"
  d$fertilizer <- stringr::str_extract(d$site, "normal|reduced")
  d$field <- dplyr::case_when(
    stringr::str_detect(d$site, "Dornburg") ~ "Dornburg",
    d$site == "Vechta" ~ "Vechta"
  )
  d <- rename_col(d, "harvest_year", "year")
  d <- rename_col(d, "grain_or_corn_dry_mass", "yield")
  d <- d[, c(1, 3:6, 8:14)]
  d <- d[d$field %in% "Dornburg", ]
  d <- d[!grepl("AFcrop", d$id), ]
  d <- rename_col(d, "id", "plot")

  d <- dplyr::bind_rows(signal16[signal16$field == "Dornburg", ], d)
  d$lat <- 51.011111
  d$lon <- 11.646111
  d$tree_species <- "poplar"
  d
}

prep_vechta <- function(file) {
  d <- read_raw(file)
  d <- d[d$site %in% "Vechta", ]
  d <- rename_col(d, "harvest_year", "year")
  d <- rename_col(d, "grain_or_corn_dry_mass", "yield")
  d$distance_to_tree_strip <- dist_from_id(d$id)
  d$data_id <- "dornburg1823"
  d$tree_species <- "poplar"
  d$yield_unit <- "g/m2"
  d <- d[, c(1, 3:6, 8:12)]
  d$field <- "Vechta"
  d <- d[!grepl("AFcrop", d$id), ]
  d <- rename_col(d, "id", "plot")
  d$lat <- 52.758056
  d$lon <- 8.534722
  d
}

prep_reiffenhausen <- function(file) {
  d <- read_raw(file)
  d <- d[, c(1:6, 11:12)]
  d <- tidyr::pivot_longer(
    d,
    cols = c(straw_2016_winter_barley, corn_2016_winter_barley,
             straw_2017_rapeseed, corn_2017_rapeseed),
    names_to = c(".value", "year", "crop"),
    names_pattern = "(straw|corn)_(\\d{4})_(.*)"
  )
  d <- rename_col(d, "corn", "yield")
  d <- rename_col(d, "straw", "yield_straw")
  d <- rename_col(d, "signal_code", "plot")
  d$distance_to_tree_strip <- dist_from_id(d$distance_from_tree_row)
  d <- drop_na_yield(d)
  d <- d[, -c(2:3)]
  d$year <- as.numeric(d$year)
  d$data_id <- "reiffenhausen16"
  d$tree_species <- "poplar"
  d$yield_unit <- "g/m2 per year"
  d$field <- "Reiffenhausen"
  d$lat <- 51.666667
  d$lon <- 9.983333
  rename_col(d, "id", "plot")
}

prep_mariensee <- function(file) {
  d <- read_raw(file, data_id = "mariensee1719")
  d <- d[, c(3:10, 15)]
  d$crop <- "grass"
  d <- rename_col(d, "grass_dm", "yield")
  d <- rename_col(d, "longitude", "lon")
  d <- rename_col(d, "latitude", "lat")
  d$tree_species <- "willow"
  d$yield_unit <- "t/ha"
  d$field <- "Mariensee"
  drop_na_yield(d)
}

prep_forst <- function(signal16, file) {
  d <- read_raw(file, data_id = "forst1920")
  d <- d[, c(2:6, 14, 19:20)]
  d <- rename_col(d, "x", "lon")
  d <- rename_col(d, "y", "lat")
  d$distance_to_tree_strip <- dist_from_id(d$signal_code)
  d <- rename_col(d, "corn_dry_weight_t_per_ha", "yield")
  d <- rename_col(d, "straw_dry_weight_t_per_ha", "yield_straw")
  d <- rename_col(d, "crop_type", "crop")
  d$year <- as.integer(format(as.Date(d$date), "%Y"))
  d$date <- NULL
  d$field <- "Forst"
  d$yield_unit <- "t/ha"
  d <- rename_col(d, "signal_code", "plot")

  sig <- signal16[signal16$field == "Forst", ]
  sig$lon <- 14.633
  sig$lat <- 51.787

  d <- dplyr::bind_rows(sig, d)
  d$tree_species <- "poplar"
  d
}

prep_gladbacherhof <- function(file) {
  d <- read_raw(file)
  d$data_id <- "gladbacherhof21"
  d$field <- "Gladbacherhof"
  d$yield_unit <- "g/m2"
  d$year <- 2021
  d <- rename_col(d, "grain_g", "yield")
  d <- rename_col(d, "biomass_g", "yield_straw")
  d <- rename_col(d, "distance", "distance_to_tree_strip")
  d <- rename_col(d, "sample", "plot")
  d$tree_species <- "mixed"
  d$crop <- "wheat"
  d
}

# Koch et al. 2025, Ihinger Hof: point measurements averaged to distance classes
prep_koch <- function(file, min_year = 2016) {
  d <- readr::read_delim(file, delim = ";", locale = readr::locale(decimal_mark = ","),
                         show_col_types = FALSE)
  lat <- mean(d$Lat)
  lon <- mean(d$Long)
  d <- janitor::clean_names(d)
  d <- rename_col(d, "treatment", "tree_species")
  d <- d |>
    dplyr::mutate(distance_to_tree_strip = dplyr::case_when(
      p_dist < 1.1  ~ 1,
      p_dist < 6.1  ~ 6,
      p_dist < 12.1  ~ 12,
      p_dist < 18.1 ~ 18,
      p_dist >= 18.1 ~ 24
    )) |>
    dplyr::group_by(year, id, block, tree_species, crop, distance_to_tree_strip) |>
    dplyr::summarise(yield = mean(yield_wweight, na.rm = TRUE), .groups = "drop")
  d$plot <- paste("Ih", d$id, d$distance_to_tree_strip, sep = "_")
  d$data_id <- "koch25"
  d$field <- "IhingerHof"
  d$lat <- lat
  d$lon <- lon
  d$yield_unit <- "t/ha"
  d[d$year >= min_year, ]
}
