# Compares the datasets of the targets pipeline with the outputs of the old
# script-based pipeline. Run after targets::tar_make():
# source(here::here("02_scriptsandmodels", "tests", "validate_against_old.R"))

library(targets)
library(dplyr)

old_project <- "C:/Users/Elizaveta/OneDrive - Universität Bayreuth/Dokumente/MasterThesis/MA_RGit_agroforestry_infield_and_beyond"
old_data <- file.path(old_project, "01_Data")

comparisons <- list(
  merged     = list(old = file.path(old_data, "AnalysisData", "AF_merged.csv"), keys = "id"),
  af_climate = list(old = file.path(old_data, "AnalysisData", "AF_climate.csv"), keys = "id"),
  af_lulc    = list(old = file.path(old_data, "AnalysisData", "AF_lulc.csv"), keys = "id"),
  af_slope   = list(old = file.path(old_data, "AnalysisData", "AF_slope.csv"), keys = "id"),
  dffinal    = list(old = tail(sort(list.files(old_data, "^dffinal_.*\\.csv$", full.names = TRUE)), 1), keys = "id"),
  moddata    = list(old = file.path(old_data, "moddata.csv"),
                    keys = c("plot", "field", "year", "crop_unified", "distance_to_tree_strip"))
)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

compare_datasets <- function(new, old, keys, name, tol = 1e-6) {
  new <- as.data.frame(new)
  old <- as.data.frame(old)
  # old moddata.csv stores the SWF matrix as SWF.swf_d100 ... or SWF.1 ... (rings of 100 m)
  names(old) <- sub("^SWF\\.(swf_d\\d+)$", "\\1", names(old))
  swf_old <- grep("^SWF\\.\\d+$", names(old), value = TRUE)
  names(old)[match(swf_old, names(old))] <- paste0("swf_d", as.numeric(sub("SWF.", "", swf_old, fixed = TRUE)) * 100)
  old <- old[!grepl("^R\\.", names(old))]

  key_new <- do.call(paste, c(lapply(new[keys], as.character), sep = "|"))
  key_old <- do.call(paste, c(lapply(old[keys], as.character), sep = "|"))
  cat("\n==", name, "==\n")
  cat("rows: new", nrow(new), "| old", nrow(old), "\n")
  cat("rows only in new:", sum(!key_new %in% key_old), "| only in old:", sum(!key_old %in% key_new), "\n")
  if (anyDuplicated(key_new) || anyDuplicated(key_old)) cat("note: keys are not unique, first match is used\n")

  only_new <- setdiff(names(new), names(old))
  only_old <- setdiff(names(old), names(new))
  if (length(only_new)) cat("columns only in new:", paste(only_new, collapse = ", "), "\n")
  if (length(only_old)) cat("columns only in old:", paste(only_old, collapse = ", "), "\n")

  m <- match(key_new, key_old)
  n <- new[!is.na(m), , drop = FALSE]
  o <- old[m[!is.na(m)], , drop = FALSE]
  cols <- setdiff(intersect(names(new), names(old)), keys)

  res <- bind_rows(lapply(cols, function(cl) {
    a <- n[[cl]]; b <- o[[cl]]
    an <- as_num(a); bn <- as_num(b)
    numeric_col <- all(is.na(an) == is.na(a)) && all(is.na(bn) == is.na(b))
    if (numeric_col) {
      d <- abs(an - bn)
      diff <- (is.na(an) != is.na(bn)) | (!is.na(d) & d > tol * pmax(1, abs(bn)))
      data.frame(dataset = name, column = cl, type = "numeric", n_compared = length(a),
                 n_different = sum(diff), max_abs_diff = suppressWarnings(max(d, na.rm = TRUE)))
    } else {
      diff <- !mapply(identical, as.character(a), as.character(b))
      data.frame(dataset = name, column = cl, type = "text", n_compared = length(a),
                 n_different = sum(diff), max_abs_diff = NA_real_)
    }
  }))
  bad <- res[res$n_different > 0, ]
  if (nrow(bad)) {
    cat("columns with differences:\n")
    print(bad[, c("column", "n_different", "max_abs_diff")], row.names = FALSE)
  } else {
    cat("all", length(cols), "common columns are identical\n")
  }
  bind_rows(res, data.frame(dataset = name, column = "<rows>", type = "rows",
                            n_compared = nrow(new), n_different = sum(!key_new %in% key_old) + sum(!key_old %in% key_new)))
}

summary_tab <- bind_rows(lapply(names(comparisons), function(nm) {
  cmp <- comparisons[[nm]]
  if (length(cmp$old) == 0 || !file.exists(cmp$old)) {
    cat("\n==", nm, "== old file not found, skipped\n")
    return(NULL)
  }
  compare_datasets(tar_read_raw(nm), readr::read_csv(cmp$old, show_col_types = FALSE, guess_max = 1e5), cmp$keys, nm)
}))

out <- file.path(res_dir <- here::here("04_resultsfigurestables", "00_validation"), "validation_against_old_pipeline.csv")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(summary_tab, out)
cat("\nDetailed comparison saved to", out, "\n")
