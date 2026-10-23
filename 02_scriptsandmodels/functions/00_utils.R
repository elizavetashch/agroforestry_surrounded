path_in <- function(...) here::here("01_inputdata", ...)

path_out <- function(...) here::here("03_outputdata", ...)

res_dir <- function(...) {
  p <- here::here("04_resultsfigurestables", ...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# accepts all spellings, e.g. "Ihinger Hof Field", "Ihinger Hof", "IhingerHof", "Ihinger"
short_field_name <- function(x) {
  s <- gsub("\\s+", "", sub("\\s*field\\s*$", "", trimws(x), ignore.case = TRUE))
  s <- ifelse(tolower(s) %in% c("ihingerhof", "ihinger"), "IhingerHof", s)
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# same formula and earth radius as geosphere::distHaversine
haversine <- function(lon1, lat1, lon2, lat2, r = 6378137) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  a <- pmin(a, 1)
  2 * atan2(sqrt(a), sqrt(1 - a)) * r
}

rename_col <- function(df, old, new) {
  names(df)[names(df) == old] <- new
  df
}

dist_from_id <- function(x) as.numeric(stringr::str_extract(x, "\\d+(?=m)"))

pick_file <- function(files, pattern) {
  f <- files[grepl(pattern, basename(files))]
  if (length(f) != 1) stop("Expected one file matching '", pattern, "', found ", length(f))
  f
}

# gives each column the type read_csv() would guess after write.csv(),
# as in the old script-based pipeline
as_csv_types <- function(df) {
  df <- as.data.frame(dplyr::ungroup(df))
  df[] <- lapply(df, as.character)
  tibble::as_tibble(readr::type_convert(df, col_types = readr::cols()))
}

save_table <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(x), file)
  file
}

save_plot <- function(p, file, width = 8, height = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = width, height = height, dpi = 300)
  file
}
