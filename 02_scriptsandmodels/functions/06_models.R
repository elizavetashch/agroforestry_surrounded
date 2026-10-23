random_terms <- "s(field, bs = \"re\") + s(crop_season, bs = \"re\") + s(year, bs = \"re\")"
dist_crop    <- "crop_unified + s(distance_to_tree_strip, k = 5, by = crop_unified)"
dist_age     <- "ti(distance_to_tree_strip, treeage, k = c(5, 5))"
swf_term     <- "s(R, by = SWF, k = 4)"
landscape    <- "l_shdi + l_contag + l_ed"

model_formulas <- function() {
  f <- function(y, ...) paste(y, "~", paste(c(...), collapse = " + "))
  data.frame(
    name = c("m_00", "m_01", "m_02", "m_03", "m_04", "m_05", "m_06", "m_07", "m_global"),
    hypothesis = c(
      "H1 distance", "H1 distance by crop", "H2 distance x tree age", "H3 SWF",
      "H4 distance + SWF", "H5 landscape", "H6 distance + SWF + landscape","H6 distance + SWF + landscape + proportion_within", "global"
    ),
    formula = c(
      f("yield_rel", "s(distance_to_tree_strip, k = 5)", random_terms),
      f("yield_rel", dist_crop, random_terms),
      f("yield_rel", dist_crop, dist_age, random_terms),
      f("yield_tha", "crop_unified", swf_term, random_terms),
      f("yield_rel", dist_crop, dist_age, swf_term, random_terms),
      f("yield_rel", "crop_unified", landscape, random_terms),
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, random_terms),
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, "prop_swf_within", random_terms),
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, 
        "treeage + AFage + PC1_c + PC1_s", random_terms)
    )
  )
}

fit_model <- function(data, spec) {
  fit <- mgcv::gam(stats::as.formula(spec$formula), data = data, method = "REML")
  list(name = spec$name, hypothesis = spec$hypothesis, formula = spec$formula, fit = fit)
}
