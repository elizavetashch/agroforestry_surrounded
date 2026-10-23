# Checks used inside the pipeline. Each returns the data unchanged, so it can be
# piped, and stops tar_make() with a clear message when a check fails.

fail <- function(name, ...) stop("[check: ", name, "] ", ..., call. = FALSE)

check_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) fail(name, "missing columns: ", paste(miss, collapse = ", "))
  invisible(df)
}

check_nrow <- function(df, n, name = "data") {
  if (nrow(df) != n) fail(name, "expected ", n, " rows, got ", nrow(df))
  invisible(df)
}

check_unique <- function(df, keys, name = "data") {
  dup <- sum(duplicated(df[keys]))
  if (dup > 0) fail(name, dup, " duplicated rows for key: ", paste(keys, collapse = " + "))
  invisible(df)
}

check_no_na <- function(df, cols, name = "data", warn = FALSE) {
  n_na <- vapply(cols, function(x) sum(is.na(df[[x]])), numeric(1))
  if (any(n_na > 0)) {
    bad <- n_na[n_na > 0]
    msg <- paste0("NA values in: ", paste0(names(bad), " (", bad, ")", collapse = ", "))
    if (warn) warning("[check: ", name, "] ", msg, call. = FALSE) else fail(name, msg)
  }
  invisible(df)
}

check_range <- function(df, col, min = -Inf, max = Inf, name = "data") {
  x <- df[[col]]
  if (any(x < min | x > max, na.rm = TRUE)) {
    fail(name, col, " outside [", min, ", ", max, "]: observed ",
         round(min(x, na.rm = TRUE), 3), " to ", round(max(x, na.rm = TRUE), 3))
  }
  invisible(df)
}

check_same_rows <- function(df, before, name = "data") {
  if (nrow(df) != nrow(before)) {
    fail(name, "row number changed from ", nrow(before), " to ", nrow(df), " (join duplicated or dropped rows)")
  }
  invisible(df)
}

field_cols <- c("yield_unit", "yield", "lat", "lon", "data_id", "field", "tree_species",
                "distance_to_tree_strip", "year", "plot", "crop")

check_field_data <- function(df, n, name) {
  df |>
    check_cols(field_cols, name) |>
    check_nrow(n, name) |>
    check_no_na(c("field", "year", "plot", "data_id", "yield_unit"), name)
}

check_moddata <- function(df, name = "moddata") {
  swf_cols <- grep("^swf_d", names(df), value = TRUE)
  df |>
    check_unique(c("plot", "field", "year", "crop_unified", "distance_to_tree_strip"), name) |>
    check_no_na(c("yield_tha", "yield_rel", "distance_to_tree_strip", "field", "year", "crop_unified"), name) |>
    check_no_na(c("treeage", "AFage", "crop_season", "PC1_c", "PC1_s", "l_shdi", "l_contag", "l_ed",
                  "prop_swf_within", swf_cols), name, warn = TRUE) |>
    check_range("yield_tha", 0, 40, name) |>
    check_range("prop_swf_within", 0, 1, name)
  for (s in swf_cols) check_range(df, s, 0, 1, name)
  m <- df |>
    dplyr::distinct(key34, key140, key140_rel) |>
    dplyr::summarise(m = mean(key140_rel), .by = key34)
  if (any(abs(m$m - 1) > 1e-8)) fail(name, "distance means do not average to 1 within field x year x crop")
  invisible(df)
}
