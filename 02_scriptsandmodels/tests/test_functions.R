# Unit tests for the helper functions, run with:
# source(here::here("02_scriptsandmodels", "tests", "test_functions.R"))
# They use small made-up data and do not need the input files.

library(dplyr)
library(ggplot2)
for (f in list.files(here::here("02_scriptsandmodels", "functions"), "\\.R$", full.names = TRUE)) source(f)

results <- c()
expect <- function(desc, cond) {
  ok <- isTRUE(tryCatch(cond, error = function(e) FALSE))
  results[desc] <<- ok
  cat(if (ok) "PASS " else "FAIL ", desc, "\n")
}
expect_error <- function(desc, expr) {
  expect(desc, inherits(tryCatch(expr, error = function(e) e), "error"))
}

# utils
expect("haversine: 1 degree latitude is about 111.3 km",
       abs(haversine(10, 50, 10, 51) - 111319.5) < 1)
expect("haversine: identical points give 0", haversine(8, 48, 8, 48) == 0)
expect("short_field_name maps long names and keeps short ones",
       identical(short_field_name(c("Ihinger Hof Field", "Forst")), c("IhingerHof", "Forst")))
expect("dist_from_id extracts the distance", identical(dist_from_id(c("D_AF_r1_24m", "D_C_r1")), c(24, NA)))
expect("as_csv_types converts character numbers", is.numeric(as_csv_types(data.frame(y = c("2016", "2017")))$y))

# merge helpers
expect("crop names are unified",
       identical(unify_crop(clean_label(c("Winter_Wheat", "sillage maize", "oat"))), c("wheat", "maize", "oat")))
expect("yield units are converted to t/ha", identical(yield_factor(c("t/ha", "g/m2", "kg/m2", "x")), c(1, 0.01, 10, NA)))
expect("SWF reference years", identical(swf_reference_year(c(2016, 2018, 2023)), c(2015, 2018, 2021)))

# checks
d <- data.frame(a = 1:3, b = c(1, NA, 3))
expect_error("check_nrow fails on wrong row number", check_nrow(d, 4))
expect_error("check_cols fails on missing column", check_cols(d, "z"))
expect_error("check_no_na fails on NA", check_no_na(d, "b"))
expect_error("check_unique fails on duplicates", check_unique(rbind(d, d), "a"))
expect_error("check_range fails outside range", check_range(d, "a", 0, 2))
expect("checks return the data when they pass", identical(check_nrow(d, 3), d))

# moddata
set.seed(1)
toy <- expand.grid(field = c("A", "B"), year = 2018:2019, crop_unified = "wheat",
                   distance_to_tree_strip = c(1, 4, 7, 24), rep = 1:3, stringsAsFactors = FALSE)
toy <- toy |>
  mutate(plot = paste(field, rep, distance_to_tree_strip, sep = "_"),
         yield_tha = runif(n(), 4, 8),
         temp_C_mean = 9 + (field == "A") + (year - 2018) / 2,
         precip_mm_sum = 600 + 50 * (field == "A") - 20 * (year - 2018),
         sun_MJ_m2_mean = 120 + 3 * (year - 2018) + (field == "B"),
         clay = ifelse(field == "A", 200, 250), sand = ifelse(field == "A", 300, 260),
         silt = ifelse(field == "A", 500, 490),
         AFage = 10, treeage = 3, l_shdi = 1, l_ed = 2, l_np = 3, l_contag = 4, mean_slope = 1,
         sowing_month = "Oct", harvest_month = "Jul") |>
  select(-rep)
swf <- data.frame(field = c("A", "B"), swf_year = 2018, prop_swf_within = 0.1,
                  swf_d100 = 0.2, swf_d200 = 0.3)
md <- make_moddata(toy, swf, c("temp_C_mean", "precip_mm_sum", "sun_MJ_m2_mean"), c("clay", "sand", "silt"))

expect("make_moddata keeps all rows", nrow(md) == nrow(toy))
expect("relative yield of the distance means averages 1 per field x year x crop",
       all(abs(tapply(md$key140_rel, md$key34, mean) - 1) < 1e-10))
expect("winter crops are classified", all(md$crop_season == "winter"))
expect("check_moddata passes on valid data", is.data.frame(suppressWarnings(check_moddata(md))))
g <- add_swf_matrix(md)
expect("SWF and R are matrices with one column per ring",
       is.matrix(g$SWF) && ncol(g$SWF) == 2 && identical(g$R[1, ], c(100, 200)))

# exploration helpers
x1 <- rnorm(100); x2 <- x1 + rnorm(100, sd = 0.1); x3 <- rnorm(100)
v <- vif(data.frame(x1, x2, x3), c("x1", "x2", "x3"))
expect("VIF detects collinear covariates", v["x1"] > 10 && v["x3"] < 2)
expect("Levene test returns a p value", levene(c(rnorm(50), rnorm(50, sd = 3)), rep(1:2, each = 50))["p_value"] < 0.05)

cat("\n", sum(results), "of", length(results), "tests passed\n")
if (!all(results)) stop("Some tests failed: ", paste(names(results)[!results], collapse = "; "))
