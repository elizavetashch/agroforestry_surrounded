# Climate ---------------------------------------------------------------------

dwd_vars <- c(
  prec_mm   = "monthly/precipitation",
  temp_C    = "monthly/air_temperature_mean",
  sun_MJ_m2 = "monthly/sun"
)

# downloads DWD monthly grids once into dir (rdwd skips files that already exist)
# and extracts the values at the field centres
extract_climate <- function(fields, years, dir) {
  gridIndex <- NULL
  utils::data("gridIndex", package = "rdwd", envir = environment())
  loc <- data.frame(x = fields$longitude, y = fields$latitude)

  out <- lapply(names(dwd_vars), function(v) {
    paths <- grep(dwd_vars[[v]], gridIndex, value = TRUE)
    ym <- stringr::str_extract(basename(paths), "\\d{6}(?=\\.asc)")
    keep <- !is.na(ym) & as.integer(substr(ym, 1, 4)) %in% years
    if (!any(keep)) stop("No DWD grids found for ", v)
    r <- rdwd::dataDWD(paths[keep], base = rdwd::gridbase, joinbf = TRUE, dir = dir, quiet = TRUE)
    r <- rdwd::projectRasterDWD(terra::rast(r), proj = "seasonal", extent = "seasonal")
    vals <- terra::extract(r, loc)[, -1, drop = FALSE]
    data.frame(
      field    = rep(fields$field, times = ncol(vals)),
      control_field = rep(fields$control_field, times = ncol(vals)),
      date     = rep(as.Date(paste0(ym[keep], "01"), "%Y%m%d"), each = nrow(vals)),
      variable = v,
      value    = unlist(vals, use.names = FALSE)
    )
  })

  dplyr::bind_rows(out) |>
    tidyr::pivot_wider(names_from = variable, values_from = value)
}

add_climate <- function(df, climate) {
  windows <- dplyr::distinct(df, field, control_field, gs_start_date, gs_end_date)
  season <- windows |>
    dplyr::left_join(climate, by = c("field", "control_field"), relationship = "many-to-many") |>
    dplyr::filter(date >= gs_start_date, date <= gs_end_date) |>
    dplyr::summarise(
      temp_C_mean    = mean(temp_C, na.rm = TRUE),
      sun_MJ_m2_mean = mean(sun_MJ_m2, na.rm = TRUE),
      precip_mm_sum  = sum(prec_mm, na.rm = TRUE),
      .by = c(field, control_field, gs_start_date, gs_end_date)
    )
  df |>
    dplyr::left_join(season, by = c("field", "control_field", "gs_start_date", "gs_end_date")) |>
    dplyr::mutate(fieldlength = haversine(min_longitude, min_latitude, max_longitude, max_latitude))
}

# Land use / land cover ---------------------------------------------------------

lulc_legend <- data.frame(
  value = c(1, 2, 4, 5, 7, 8, 9, 10, 11),
  class = c("Water", "Trees", "Flooded Vegetation", "Crops", "Built Area",
            "Bare Ground", "Snow/Ice", "Clouds", "Rangeland")
)

# landscape metrics within a 3 km buffer around each field centre;
# Forst lies in UTM zone 33U, all other fields in 32U
lulc_metrics <- function(fieldpolygons_file, lulc_files, years = 2017:2023, buffer = 3000) {
  fp <- read_fieldpolygons(fieldpolygons_file)
  out <- list()
  for (zone in c("32U", "33U")) {
    sites <- if (zone == "33U") fp[fp$field == "Forst", ] else fp[fp$field != "Forst", ]
    for (yr in years) {
      f <- lulc_files[basename(dirname(lulc_files)) == yr & startsWith(basename(lulc_files), paste0(zone, "_"))]
      if (!length(f)) stop("No LULC raster for zone ", zone, " year ", yr)
      r <- terra::rast(f[1])
      levels(r) <- lulc_legend
      buf <- sf::st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326) |>
        sf::st_transform(terra::crs(r)) |>
        sf::st_buffer(dist = buffer)
      for (i in seq_len(nrow(buf))) {
        p <- terra::vect(buf[i, ])
        rc <- terra::mask(terra::crop(r, p), p)
        out[[length(out) + 1]] <- data.frame(
          field    = sites$field[i],
          control_field = sites$control_field[i],
          year     = yr,
          l_contag = landscapemetrics::lsm_l_contag(rc)$value,
          l_np     = landscapemetrics::lsm_l_np(rc)$value,
          l_sidi   = landscapemetrics::lsm_l_sidi(rc)$value,
          l_shdi   = landscapemetrics::lsm_l_shdi(rc)$value,
          l_ai     = landscapemetrics::lsm_l_ai(rc)$value,
          l_ed     = landscapemetrics::lsm_l_ed(rc)$value
        )
      }
    }
  }
  res <- dplyr::bind_rows(out)
  # no satellite coverage in 2016: 2017 values are used
  res16 <- dplyr::mutate(dplyr::filter(res, year == min(years)), year = 2016)
  dplyr::arrange(dplyr::bind_rows(res, res16), field, control_field, year)
}

add_lulc <- function(df, lulc) {
  dplyr::left_join(df, lulc, by = c("field", "control_field", "year"), relationship = "many-to-one")
}

# Slope (mean and SD of the DTM slope within each field polygon, columns MEAN and STD
# of fieldpolygons_with_control.csv) ---------------------------------------------

add_slope <- function(df, fieldpolygons_file) {
  slope <- read_fieldpolygons(fieldpolygons_file) |>
    dplyr::select(field, control_field, mean_slope, sd_slope)
  dplyr::left_join(df, slope, by = c("field", "control_field"), relationship = "many-to-one")
}

# Soil texture (SoilGrids, sampled at the field centre) -------------------------

soil_texture <- function(fieldpolygons_file, soil_files, soil_vars = c("clay", "sand", "silt")) {
  r <- terra::rast(soil_files)
  for (v in soil_vars) {
    hit <- grepl(v, names(r), ignore.case = TRUE)
    if (!v %in% names(r) && sum(hit) == 1) names(r)[hit] <- v
  }
  if (!all(soil_vars %in% names(r))) {
    stop("SoilGrids layers are named: ", paste(names(r), collapse = ", "),
         ". Expected one layer each for: ", paste(soil_vars, collapse = ", "))
  }
  fp <- read_fieldpolygons(fieldpolygons_file)
  pts <- sf::st_as_sf(fp, coords = c("Longitude", "Latitude"), crs = 4326) |>
    sf::st_transform(terra::crs(r))
  vals <- terra::extract(r, sf::st_coordinates(pts), method = "simple")
  vals$ID <- NULL
  dplyr::bind_cols(data.frame(field = fp$field, control_field = fp$control_field), vals)
}

# Small woody features (SWF) ------------------------------------------------------

# raw SWF maps per reference year, read from 01_inputdata/swf_raw/<year>/:
#   2015: swf_2015_mosaic.tif, class 3 = woody feature (made binary here)
#   2018, 2021: original tiles (0/1), combined into one virtual raster
swf_raw_raster <- function(swf_raw_files, yr) {
  f <- swf_raw_files[basename(dirname(swf_raw_files)) == yr]
  if (!length(f)) stop("No raw SWF raster for ", yr)
  if (length(f) == 1) terra::rast(f) else terra::vrt(f)
}

swf_binary <- function(r, yr) {
  if (yr == "2015") terra::ifel(r == 3, 1, 0) else as.numeric(r)
}

swf_share <- function(raster, geom) {
  vals <- terra::extract(raster, terra::vect(geom))[[2]]
  vals <- vals[!is.na(vals)]
  if (length(vals)) sum(vals) / length(vals) else NA_real_
}

# share of woody pixels in rings of 100 m width outside the field boundary,
# and inside the field (disk with half the field diagonal as radius),
# for agroforestry and control fields; each map is cut to a 3 km buffer first
swf_profiles <- function(fieldpolygons_file, swf_raw_files, distances = seq(100, 1000, by = 100),
                         years = c("2015", "2018", "2021"), buffer = 3000) {
  fp <- read_fieldpolygons(fieldpolygons_file)
  fp$boundary_min <- haversine(fp$minLongitude, fp$minLatitude, fp$maxLongitude, fp$maxLatitude) / 2

  out <- list()
  for (yr in years) {
    r_all <- swf_raw_raster(swf_raw_files, yr)
    pts <- sf::st_transform(sf::st_as_sf(fp, coords = c("Longitude", "Latitude"), crs = 4326), terra::crs(r_all))
    for (i in seq_len(nrow(pts))) {
      g <- sf::st_geometry(pts[i, ])
      buf <- terra::vect(sf::st_buffer(g, buffer))
      r <- swf_binary(terra::mask(terra::crop(r_all, buf), buf), yr)
      bmin <- pts$boundary_min[i]
      rings <- vapply(seq_along(distances), function(k) {
        inner <- if (k == 1) bmin else bmin + distances[k - 1]
        swf_share(r, sf::st_difference(sf::st_buffer(g, bmin + distances[k]), sf::st_buffer(g, inner)))
      }, numeric(1))
      out[[length(out) + 1]] <- data.frame(
        field = pts$field[i],
        control_field = pts$control_field[i],
        swf_year = as.numeric(yr),
        prop_swf_within = swf_share(r, sf::st_buffer(g, bmin)),
        t(stats::setNames(rings, paste0("swf_d", distances)))
      )
    }
  }
  dplyr::arrange(dplyr::bind_rows(out), field, control_field, swf_year)
}
