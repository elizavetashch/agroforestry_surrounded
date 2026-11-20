random_terms <- "s(key34, bs = \"re\")"

model_terms_library <- c(
  crop      = "crop_unified",
  site      = "treeage + AFage + PC1_c + PC1_s + sd_slope",
  site_ctrl = "PC1_c + PC1_s + sd_slope",
  beyond_param = "s(swf_beyond, k = 3) + l_shdi + l_contag + l_ed",
  beyond    = "s(swf_beyond, k = 3) + s(l_shdi, k = 3) + l_contag + l_ed",
  beyond_fn = "s(R, by = SWF, k = 4) + s(l_shdi, k = 3) + l_contag + l_ed",
  dist      = "s(distance_to_tree_strip, k = 5)",
  dist_crop = "s(distance_to_tree_strip, k = 5, by = crop_unified)",
  dist_age  = "ti(distance_to_tree_strip, treeage, k = c(5, 5))",
  infield   = "s(prop_swf_within, k = 3)",
  int_swf   = "ti(swf_beyond, prop_swf_within, k = c(3, 3))",
  int_shdi  = "ti(l_shdi, prop_swf_within, k = c(3, 3))"
)

# control = TRUE: same models for the control fields, without the terms that need trees
# (distance, distance x tree age, tree age, AF age); models that become identical are dropped
model_formulas <- function(response = "yield_tha", crop_main = TRUE, control = FALSE) {
  f <- function(keys) {
    if (!crop_main) keys <- setdiff(keys, "crop")
    if (control) {
      keys <- setdiff(keys, c("dist", "dist_crop", "dist_age"))
      keys[keys == "site"] <- "site_ctrl"
    }
    paste(response, "~", paste(c(model_terms_library[keys], random_terms), collapse = " + "))
  }
  h1 <- c("crop", "site", "beyond")
  h1p <- c("crop", "site", "beyond_param")
  h2 <- c(h1, "dist_crop", "dist_age", "infield")
  h2p <- c(h1p, "dist_crop", "dist_age", "infield")
  
  specs <- list(
    m00  = list(h = "H1 baseline: regional model without in-field heterogeneity", keys = h1),
    m00p  = list(h = "H1 baseline: regional model without in-field heterogeneity", keys = h1p),
    m01  = list(h = "H1: + distance to tree strip", keys = c(h1, "dist")),
    m01p  = list(h = "H1: + distance to tree strip", keys = c(h1p, "dist")),
    m01b = list(h = "H1: + distance to tree strip by crop", keys = c(h1, "dist_crop")),
    m01bp = list(h = "H1: + distance to tree strip by crop", keys = c(h1p, "dist_crop")),
    m01c = list(h = "H1: + distance by crop + distance x tree age", keys = c(h1, "dist_crop", "dist_age")),
    m01cp = list(h = "H1: + distance by crop + distance x tree age", keys = c(h1p, "dist_crop", "dist_age")),
    m01d = list(h = "SWF: m01c with the ring-resolved (functional) SWF effect instead of mean SWF",
                keys = c("crop", "site", "beyond_fn", "dist_crop", "dist_age")),
    m01dp = list(h = "SWF: m01c with the ring-resolved (functional) SWF effect instead of mean SWF",
                keys = c("crop", "site", "beyond_param", "dist_crop", "dist_age")),
    m02  = list(h = "H2 baseline: m01c + in-field complexity, no interaction", keys = h2),
    m02p  = list(h = "H2 baseline: m01c + in-field complexity, no interaction", keys = h2p),
    m02b = list(h = "H2: SWF beyond x in-field complexity", keys = c(h2, "int_swf")),
    m02bp = list(h = "H2: SWF beyond x in-field complexity", keys = c(h2p, "int_swf")),
    m02c = list(h = "H2: Shannon diversity x in-field complexity", keys = c(h2, "int_shdi")),
    m02cp = list(h = "H2: Shannon diversity x in-field complexity", keys = c(h2p, "int_shdi")),
    m02d = list(h = "H2: SWF beyond and Shannon diversity x in-field complexity", keys = c(h2, "int_swf", "int_shdi")),
    m02dp = list(h = "H2: SWF beyond and Shannon diversity x in-field complexity", keys = c(h2p, "int_swf", "int_shdi"))
  )
  
  out <- data.frame(
    name       = names(specs),
    hypothesis = vapply(specs, `[[`, character(1), "h"),
    formula    = vapply(specs, function(s) f(s$keys), character(1)),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  if (control) out[!duplicated(out$formula), ] else out
}

fit_model <- function(data, spec) {
  f <- stats::as.formula(spec$formula)
  list(
    name       = spec$name,
    hypothesis = spec$hypothesis,
    formula    = spec$formula,
    fit_ml     = mgcv::gam(f, data = data, family = Gamma(link = "log"), method = "ML"),
    fit        = mgcv::gam(f, data = data, family = Gamma(link = "log"), method = "REML")
  )
}