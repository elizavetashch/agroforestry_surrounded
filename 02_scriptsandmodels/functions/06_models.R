library(mgcv)

random_terms <- "s(field, bs = \"re\") + s(year, bs = \"re\") + s(key34, bs = \"re\")"

model_formulas <- function() {
  f <- function(y, ...) paste(y, "~", paste(c(...), collapse = " + "))
  
  # --- Fixed-effect building blocks -----------------------------------------
  terms <- c(
    crop      = "crop_unified",
    dist      = "s(distance_to_tree_strip, k = 5)",
    dist_crop = "crop_unified + s(distance_to_tree_strip, k = 5, by = crop_unified)",
    dist_age  = "ti(distance_to_tree_strip, treeage, k = c(5, 5))",
    swf       = "s(R, by = SWF, k = 4)",
    land      = "l_shdi + l_contag + l_ed",
    prop_in   = "prop_swf_within",
    site      = "treeage + AFage + PC1_c + PC1_s"
  )
  
  # --- Model specifications --------------------------------------------------
  specs <- list(
    m_00     = list(h = "H1a: yield varies with distance to the tree strip",
                    keys = c("crop", "dist")),
    m_01     = list(h = "H1b: the distance effect differs among crops",
                    keys = c("dist_crop")),
    m_02     = list(h = "H2: the distance effect changes with tree age",
                    keys = c("dist_crop", "dist_age")),
    m_03     = list(h = "H3: surrounding SWF affect yield (no distance effect)",
                    keys = c("crop", "swf")),
    m_04     = list(h = "H4: surrounding SWF improve the distance model",
                    keys = c("dist_crop", "dist_age", "swf")),
    m_05     = list(h = "H5: landscape composition and configuration affect yield (no distance effect)",
                    keys = c("crop", "land")),
    m_06     = list(h = "H6: landscape improves the distance + SWF model",
                    keys = c("dist_crop", "dist_age", "swf", "land")),
    m_07     = list(h = "H7: in-field SWF proportion adds information beyond landscape and surrounding SWF",
                    keys = c("dist_crop", "dist_age", "swf", "land", "prop_in")),
    m_08     = list(h = "H8: in-field SWF proportion alone improves the distance model",
                    keys = c("dist_crop", "dist_age", "prop_in")),
    m_09     = list(h = "H9: surrounding SWF are a better spill-over proxy than landscape heterogeneity",
                    keys = c("dist_crop", "dist_age", "swf", "prop_in")),
    m_10     = list(h = "H10: site characteristics improve the model",
                    keys = c("dist_crop", "dist_age", "swf", "land", "site")),
    m_global = list(h = "Global: in-field complexity effect shaped by landscape, controlling for site",
                    keys = c("dist_crop", "dist_age", "swf", "land", "site", "prop_in"))
  )
  
  data.frame(
    name       = names(specs),
    hypothesis = vapply(specs, `[[`, character(1), "h"),
    formula    = vapply(specs, function(s) f("yield_tha", terms[s$keys], random_terms),
                        character(1)),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
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
