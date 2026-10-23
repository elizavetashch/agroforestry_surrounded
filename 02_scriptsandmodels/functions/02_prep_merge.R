crop_lookup <- c(
  "sillage maize"      = "maize",
  "silage maize"       = "maize",
  "winterwheat"        = "wheat",
  "ww"                 = "wheat",
  "winter wheat"       = "wheat",
  "wb"                 = "barley",
  "wp"                 = "pea",
  "winter barley"      = "barley",
  "spring barley"      = "barley",
  "summer barley"      = "barley",
  "wr"                 = "rapeseed",
  "winter rye"         = "rye",
  "oil rape"           = "rapeseed",
  "winter oilseedrape" = "rapeseed",
  "winter rapeseed"    = "rapeseed",
  "bristle oat"        = "oat"
)

clean_label <- function(x) {
  x |> stringr::str_replace_all("[-_]", " ") |> stringr::str_to_lower() |> stringr::str_squish()
}

unify_crop <- function(crop) {
  dplyr::if_else(crop %in% names(crop_lookup), unname(crop_lookup[crop]), crop)
}

yield_factor <- function(unit) {
  dplyr::case_when(
    unit == "t/ha"                        ~ 1,
    unit == "kg/m2"                       ~ 10,
    unit %in% c("g/m2", "g/m2 per year")  ~ 0.01,
    .default = NA_real_
  )
}

merge_fields <- function(fields, fieldpolygons_file, fertilization_file) {
  d <- dplyr::bind_rows(lapply(fields, function(x) {
    dplyr::mutate(as_csv_types(x), dplyr::across(dplyr::any_of(c("plot", "crop", "tree_species", "fertilizer")), as.character))
  })) |>
    dplyr::mutate(
      crop            = clean_label(crop),
      tree_species    = clean_label(tree_species),
      crop_unified    = unify_crop(crop),
      yield_tha       = yield * yield_factor(yield_unit),
      yield_straw_tha = yield_straw * yield_factor(yield_unit),
      harvestyear = dplyr::case_when(
        field == "Dornburg"                  ~ 2014,
        field == "Forst"                     ~ 2014,
        field == "Gladbacherhof"             ~ 2019,
        field == "IhingerHof" & year < 2019  ~ 2016,
        field == "IhingerHof" & year >= 2022 ~ 2022,
        field == "IhingerHof" & year >= 2019 ~ 2019,
        field == "Mariensee"                 ~ 2015,
        field == "Reiffenhausen"             ~ 2015,
        field == "Vechta"                    ~ 2018,
        field == "Wendhausen" & year < 2021  ~ 2013,
        field == "Wendhausen" & year >= 2021 ~ 2021
      ),
      treeage = year - harvestyear
    ) |>
    dplyr::filter(!is.na(distance_to_tree_strip), !is.na(yield)) |>
    dplyr::filter(!(year == 2020 & plot == "V_AF_r4_4m"),
                  !(year == 2021 & plot == "V_AF_r1_18m")) |>
    dplyr::mutate(plot = dplyr::if_else(plot == "F_AF_r3_24m" & dplyr::near(yield, 625.42),
                                        "F_AF_r4_24m", plot)) |>
    dplyr::filter(!is.na(yield_tha))

  polygons <- readr::read_csv(fieldpolygons_file, show_col_types = FALSE) |>
    dplyr::mutate(field = short_field_name(Name)) |>
    dplyr::select(field, Latitude, Longitude, minLatitude, minLongitude,
                  maxLatitude, maxLongitude, Area) |>
    dplyr::distinct()

  d <- d |>
    dplyr::left_join(polygons, by = "field") |>
    dplyr::mutate(id = dplyr::case_when(
      stringr::str_detect(plot, "m$")     ~ paste0(plot, "_", year),
      stringr::str_detect(plot, "^Ih")    ~ paste0(plot, "m_", year),
      stringr::str_detect(plot, "^[A-Z]") ~ paste0(plot, "_", distance_to_tree_strip, "m_", year),
      field == "Gladbacherhof"            ~ paste0(plot, "_", distance_to_tree_strip, "m_", year)
    ))

  dup <- duplicated(d$id) | duplicated(d$id, fromLast = TRUE)
  if (any(dup)) message("merge_fields: removing ", sum(dup), " rows with duplicated id")
  d <- d[!dup, ]

  fertilization <- readr::read_csv(fertilization_file, show_col_types = FALSE) |>
    janitor::clean_names()

  d$fertilizer[is.na(d$fertilizer)] <- "normal"

  d |>
    dplyr::left_join(fertilization, by = c("field", "year" = "harvest_year", "crop_unified", "fertilizer")) |>
    dplyr::mutate(
      growing_start = match(sowing_month, month.abb),
      growing_end   = match(harvest_month, month.abb),
      gs_start_year = dplyr::if_else(growing_start >= 7, year - 1, year),
      gs_start_date = as.Date(sprintf("%04d-%02d-01", as.integer(gs_start_year), growing_start)),
      gs_end_date   = as.Date(sprintf("%04d-%02d-01", as.integer(year), growing_end))
    ) |>
    dplyr::select(-gs_start_year) |>
    janitor::clean_names()
}
