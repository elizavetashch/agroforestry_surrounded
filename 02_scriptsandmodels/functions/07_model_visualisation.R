model_comparison <- function(models) {
  dplyr::bind_rows(lapply(models, function(m) {
    s <- summary(m$fit)
    data.frame(
      model = m$name, hypothesis = m$hypothesis, response = all.vars(m$fit$formula)[1],
      n = s$n, df = attr(stats::logLik(m$fit), "df"), AIC = stats::AIC(m$fit),
      r2_adj = s$r.sq, deviance_explained = s$dev.expl, formula = m$formula
    )
  })) |>
    dplyr::mutate(delta_AIC = AIC - min(AIC), .after = AIC, .by = response) |>
    dplyr::arrange(response, AIC)
}

model_terms <- function(models) {
  dplyr::bind_rows(lapply(models, function(m) {
    s <- summary(m$fit)
    p <- data.frame(model = m$name, type = "parametric", term = rownames(s$p.table),
                    estimate = s$p.table[, 1], se = s$p.table[, 2], statistic = s$p.table[, 3],
                    p_value = s$p.table[, 4], row.names = NULL)
    sm <- data.frame(model = m$name, type = "smooth", term = rownames(s$s.table),
                     edf = s$s.table[, "edf"], ref_df = s$s.table[, "Ref.df"],
                     statistic = s$s.table[, 3], p_value = s$s.table[, 4], row.names = NULL)
    dplyr::bind_rows(p, sm)
  }))
}

model_kcheck <- function(models) {
  dplyr::bind_rows(lapply(models, function(m) {
    k <- mgcv::k.check(m$fit)
    data.frame(model = m$name, term = rownames(k), k, row.names = NULL, check.names = FALSE)
  }))
}

save_model_tables <- function(models, dir) {
  c(save_table(model_comparison(models), file.path(dir, "model_comparison_aic.csv")),
    save_table(model_terms(models), file.path(dir, "model_terms.csv")),
    save_table(model_kcheck(models), file.path(dir, "model_kcheck.csv")))
}

save_model_figures <- function(models, dir) {
  unlist(lapply(models, function(m) {
    f_check <- file.path(dir, paste0(m$name, "_diagnostics.png"))
    grDevices::png(f_check, width = 2400, height = 2400, res = 300)
    graphics::par(mfrow = c(2, 2))
    utils::capture.output(mgcv::gam.check(m$fit, pch = 16, cex = 0.5))
    grDevices::dev.off()

    f_smooth <- file.path(dir, paste0(m$name, "_smooths.png"))
    grDevices::png(f_smooth, width = 3000, height = 3000, res = 300)
    plot(m$fit, pages = 1, shade = TRUE, shade.col = "lightblue", seWithMean = TRUE,
         scale = 0, residuals = TRUE, pch = 16, cex = 0.3, col = "grey50")
    grDevices::dev.off()

    c(f_check, f_smooth)
  }))
}

plot_relative_yield <- function(moddata, file) {
  means <- moddata |>
    dplyr::summarise(yield_tha = mean(yield_tha, na.rm = TRUE),
                     .by = c(field, year, crop_unified, distance_to_tree_strip))
  p <- ggplot(moddata, aes(distance_to_tree_strip, yield_tha)) +
    geom_jitter(aes(colour = field), width = 0.3, height = 0, alpha = 0.25, size = 1.2) +
    geom_line(data = means, aes(group = interaction(field, year), colour = field),
              alpha = 0.6, linewidth = 0.4) +
    stat_summary(fun.data = mean_se, geom = "pointrange", colour = "black", size = 0.3) +
    facet_wrap(~ crop_unified, scales = "free_y") +
    labs(x = "Distance to tree strip (m)", y = expression(Yield~(t~ha^-1)), colour = "Site") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  save_plot(p, file, 10, 7)
}

plot_distance_by_crop <- function(models, data, file, model = "m_01") {
  m <- models[[which(vapply(models, `[[`, "", "name") == model)]]$fit
  crops <- levels(data$crop_unified)
  nd <- expand.grid(distance_to_tree_strip = seq(min(data$distance_to_tree_strip),
                                                 max(data$distance_to_tree_strip), length.out = 100),
                    crop_unified = crops)
  for (v in intersect(c("key34", "field", "year", "crop_season"), names(data))) {
    nd[[v]] <- data[[v]][!is.na(data[[v]])][1]
  }
  ex <- grep("key34|field|year|crop_season", vapply(m$smooth, `[[`, "", "label"), value = TRUE)
  pr <- stats::predict(m, nd, type = "link", se.fit = TRUE, exclude = ex)
  nd$fit <- exp(pr$fit)
  nd$lwr <- exp(pr$fit - 1.96 * pr$se.fit)
  nd$upr <- exp(pr$fit + 1.96 * pr$se.fit)
  rng <- dplyr::summarise(data, lo = min(distance_to_tree_strip), hi = max(distance_to_tree_strip), .by = crop_unified)
  nd <- dplyr::left_join(nd, rng, by = "crop_unified") |>
    dplyr::filter(distance_to_tree_strip >= lo, distance_to_tree_strip <= hi)
  p <- ggplot(nd, aes(distance_to_tree_strip, fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80") +
    geom_line() +
    facet_wrap(~ crop_unified, scales = "free_y") +
    labs(x = "Distance to tree strip (m)", y = bquote("Predicted yield (t ha"^-1*"), " * .(model))) +
    theme_bw()
  save_plot(p, file, 9, 6)
}

