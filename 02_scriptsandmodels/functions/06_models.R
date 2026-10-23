random_terms <- "s(field, bs = \"re\") + s(crop_season, bs = \"re\") + s(year, bs = \"re\")"
dist_crop    <- "crop_unified + s(distance_to_tree_strip, k = 5, by = crop_unified)"
dist_age     <- "ti(distance_to_tree_strip, treeage, k = c(5, 5))"
swf_term     <- "s(R, by = SWF, k = 4)"
landscape    <- "l_shdi + l_contag + l_ed"

model_formulas <- function() {
  f <- function(y, ...) paste(y, "~", paste(c(...), collapse = " + "))
  data.frame(
    name = c("m_00", "m_01", "m_02", "m_03", "m_04", "m_05", "m_06", "m_07", "m_08", "m_09", "m_10", "m_global"),
    hypothesis = c(
      "H1 distance", "H1 distance by crop", "H2 distance by crop x tree age", 
      "H3 SWF + no distance, but with crop effect",
      "H4: H2 + SWF", "H5 landscape + no distance, but with crop effect", 
      "H6  distance by  + dist_age + SWF + landscape",
      "H7 distance by  + dist_age + SWF + landscape + proportion_within", 
      "H8: distance by  + dist_age + proportion_within", 
      "H9: distance by  + dist_age + SWF + landscape + proportion_within", 
      "H10: dist_crop, dist_age, swf_term, landscape, treeage + AFage + PC1_c + PC1_s",
      "global"
    ),
    formula = c(
      f("yield_rel", "s(distance_to_tree_strip, k = 5)", random_terms), #m_00
      f("yield_rel", dist_crop, random_terms), #m_01
      f("yield_rel", dist_crop, dist_age, random_terms), #m_02
      f("yield_tha", "crop_unified", swf_term, random_terms), #m_03
      f("yield_rel", dist_crop, dist_age, swf_term, random_terms), #m_04
      f("yield_rel", "crop_unified", landscape, random_terms), #m_05
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, random_terms), #m_06
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, "prop_swf_within", random_terms), #m_07
      f("yield_rel", dist_crop, dist_age, "prop_swf_within", random_terms), #m_08
      f("yield_rel", dist_crop, dist_age, swf_term, "prop_swf_within", random_terms), #m_09
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, #m_10
        "treeage + AFage + PC1_c + PC1_s", random_terms),
      f("yield_rel", dist_crop, dist_age, swf_term, landscape, #m_global
        "treeage + AFage + PC1_c + PC1_s", random_terms, "prop_swf_within")
    )
  )
}

fit_model <- function(data, spec) {
  fit <- mgcv::gam(stats::as.formula(spec$formula), data = data, method = "REML")
  list(name = spec$name, hypothesis = spec$hypothesis, formula = spec$formula, fit = fit)
}
