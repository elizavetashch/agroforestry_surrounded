model_comparison <- function(models) {
  dplyr::bind_rows(lapply(models, function(m) {
    s <- summary(m$fit)
    data.frame(
      model = m$name, hypothesis = m$hypothesis, response = all.vars(m$fit$formula)[1],
      n = s$n,
      df_ML = attr(stats::logLik(m$fit_ml), "df"),
      AIC_ML = stats::AIC(m$fit_ml),
      r2_adj_REML = s$r.sq,
      deviance_explained_REML = s$dev.expl,
      formula = m$formula
    )
  })) |>
    dplyr::mutate(delta_AIC = AIC_ML - min(AIC_ML), .after = AIC_ML, .by = response) |>
    dplyr::arrange(response, AIC_ML)
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


# readable names for model terms -----------------------------------------------

term_labels <- c(
  "treeage" = "Tree age", "AFage" = "Agroforestry system age", "PC1_c" = "PC1 climate",
  "PC1_s" = "PC1 soil", "l_contag" = "Contagion", "l_ed" = "Edge density",
  "s(swf_beyond)" = "SWF mean", "s(l_shdi)" = "Shannon diversity",
  "s(prop_swf_within)" = "In-field complexity", "s(R):SWF" = "SWF rings (functional)",
  "s(distance_to_tree_strip)" = "Distance", "ti(distance_to_tree_strip,treeage)" = "Distance x tree age",
  "ti(swf_beyond,prop_swf_within)" = "SWF mean x in-field complexity",
  "ti(l_shdi,prop_swf_within)" = "Shannon x in-field complexity", "s(key34)" = "Random: field x year x crop"
)

nice_term <- function(x) {
  out <- ifelse(x %in% names(term_labels), term_labels[x], x)
  out <- sub("^crop_unified", "Crop: ", out)
  out <- sub("^s\\(distance_to_tree_strip\\):crop_unified", "Distance: ", out)
  unname(out)
}

sig_class <- function(p) {
  cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("p < 0.001", "p < 0.01", "p < 0.05", "n.s."))
}

# parametric effects of all models as % change in yield with 95% CI -----------
plot_coef_forest <- function(models, file) {
  lv <- vapply(models, `[[`, "", "name")
  d <- model_terms(models) |>
    dplyr::filter(type == "parametric", term != "(Intercept)") |>
    dplyr::mutate(
      effect = 100 * (exp(estimate) - 1),
      lwr = 100 * (exp(estimate - 1.96 * se) - 1),
      upr = 100 * (exp(estimate + 1.96 * se) - 1),
      significant = ifelse(p_value < 0.05, "p < 0.05", "n.s."),
      term = nice_term(term),
      model = factor(model, levels = rev(lv))
    )
  p <- ggplot(d, aes(effect, model, colour = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_point(size = 2) +
    facet_wrap(~ term, scales = "free_x") +
    scale_colour_manual(values = c("p < 0.05" = "black", "n.s." = "grey65")) +
    labs(x = "Effect on yield (%) with 95% CI", y = NULL, colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  save_plot(p, file, 12, 9)
}

# smooth terms of all models: significance and edf ---------------------------
plot_smooth_heatmap <- function(models, file) {
  lv <- vapply(models, `[[`, "", "name")
  d <- model_terms(models) |>
    dplyr::filter(type == "smooth") |>
    dplyr::mutate(term = nice_term(term), significance = sig_class(p_value),
                  model = factor(model, levels = lv))
  order_terms <- unique(d$term)
  order_terms <- rev(c(setdiff(order_terms, "Random: field x year x crop"), "Random: field x year x crop"))
  d$term <- factor(d$term, levels = order_terms)
  p <- ggplot(d, aes(model, term, fill = significance)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.1f", edf)), size = 3) +
    scale_fill_manual(values = c("p < 0.001" = "grey25", "p < 0.01" = "grey50",
                                 "p < 0.05" = "grey75", "n.s." = "grey95"), drop = FALSE) +
    labs(x = NULL, y = NULL, fill = NULL, caption = "Numbers = effective degrees of freedom (1 = linear)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid = element_blank())
  save_plot(p, file, 9, 7)
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

# typical-value row: numeric covariates at their median, other columns from the first row
typical_rows <- function(data, n) {
  nd <- data[rep(1, n), ]
  for (v in names(data)) {
    if (is.numeric(data[[v]]) && !is.matrix(data[[v]])) nd[[v]] <- stats::median(data[[v]], na.rm = TRUE)
  }
  nd
}

random_labels <- function(m) {
  grep("key34|field|year|crop_season", vapply(m$smooth, `[[`, "", "label"), value = TRUE)
}

plot_distance_by_crop <- function(models, data, file, model = "m01c") {
  m <- models[[which(vapply(models, `[[`, "", "name") == model)]]$fit
  grid <- expand.grid(distance_to_tree_strip = seq(min(data$distance_to_tree_strip),
                                                   max(data$distance_to_tree_strip), length.out = 100),
                      crop_unified = levels(data$crop_unified))
  nd <- typical_rows(data, nrow(grid))
  nd$distance_to_tree_strip <- grid$distance_to_tree_strip
  nd$crop_unified <- factor(grid$crop_unified, levels = levels(data$crop_unified))
  pr <- stats::predict(m, nd, type = "link", se.fit = TRUE, exclude = random_labels(m))
  grid$fit <- exp(pr$fit)
  grid$lwr <- exp(pr$fit - 1.96 * pr$se.fit)
  grid$upr <- exp(pr$fit + 1.96 * pr$se.fit)
  rng <- dplyr::summarise(data, lo = min(distance_to_tree_strip), hi = max(distance_to_tree_strip), .by = crop_unified)
  grid <- dplyr::left_join(grid, rng, by = "crop_unified") |>
    dplyr::filter(distance_to_tree_strip >= lo, distance_to_tree_strip <= hi)
  p <- ggplot(grid, aes(distance_to_tree_strip, fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80") +
    geom_line() +
    facet_wrap(~ crop_unified, scales = "free_y") +
    labs(x = "Distance to tree strip (m)",
         y = bquote("Predicted yield (t ha"^-1*"), " * .(model) * ", other covariates at median")) +
    theme_bw()
  save_plot(p, file, 9, 6)
}

plot_aic <- function(models, file) {
  d <- model_comparison(models) |>
    dplyr::mutate(model = factor(model, levels = rev(vapply(models, `[[`, "", "name"))),
                  hypothesis_group = dplyr::case_when(
                    model == "m01d"              ~ "SWF functional effect",
                    grepl("^m00|^m01", model)    ~ "H1: in-field heterogeneity",
                    .default                     = "H2: beyond x in-field interaction"))
  p <- ggplot(d, aes(delta_AIC, model, fill = hypothesis_group)) +
    geom_col() +
    geom_vline(xintercept = 2, linetype = "dashed", colour = "grey40") +
    labs(x = expression(Delta * "AIC (ML) relative to the best model"), y = NULL, fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  save_plot(p, file, 7, 4.5)
}

# functional SWF effect beta(r): weight of woody features in each 100 m ring
plot_swf_function <- function(models, data, file, model = "m01d") {
  m <- models[[which(vapply(models, `[[`, "", "name") == model)]]$fit
  r <- seq(min(data$R), max(data$R), length.out = 100)
  nd <- typical_rows(data, length(r))
  nd$R <- matrix(r, ncol = 1)
  nd$SWF <- matrix(1, nrow = length(r), ncol = 1)
  X <- stats::predict(m, nd, type = "lpmatrix")
  sm <- m$smooth[[which(vapply(m$smooth, `[[`, "", "label") == "s(R):SWF")]]
  cols <- sm$first.para:sm$last.para
  d <- data.frame(r = r, beta = drop(X[, cols] %*% stats::coef(m)[cols]),
                  se = sqrt(rowSums((X[, cols] %*% m$Vp[cols, cols]) * X[, cols])))
  p <- ggplot(d, aes(r, beta)) +
    geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), fill = "grey80") +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(x = "Distance of ring from field boundary (m)",
         y = expression(beta(r) ~ "(effect of woody share on log yield)")) +
    theme_bw()
  save_plot(p, file, 7, 4.5)
}

# H2 -------------------------------------------------------------------------

plot_heterogeneity_space <- function(moddata, beyond, infield, file) {
  d <- dplyr::distinct(moddata, field, year, .data[[beyond]], .data[[infield]])
  p <- ggplot(d, aes(.data[[beyond]], .data[[infield]], colour = field)) +
    geom_point(size = 3) +
    labs(x = paste("Beyond-field heterogeneity:", beyond),
         y = paste("In-field complexity:", infield), colour = "Site") +
    theme_bw()
  save_plot(p, file, 7, 5)
}

plot_interaction_raw <- function(moddata, beyond, infield, file) {
  d <- moddata |>
    dplyr::mutate(yield_crop = yield_tha / mean(yield_tha), .by = crop_unified) |>
    dplyr::summarise(yield_crop = mean(yield_crop),
                     .by = c(field, year, crop_unified, dplyr::all_of(c(beyond, infield)))) |>
    dplyr::mutate(infield_class = ifelse(.data[[infield]] > stats::median(.data[[infield]]),
                                         "high in-field complexity", "low in-field complexity"))
  p <- ggplot(d, aes(.data[[beyond]], yield_crop, colour = infield_class)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
    labs(x = paste("Beyond-field heterogeneity:", beyond),
         y = "Yield relative to crop mean (field x year x crop means)", colour = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  save_plot(p, file, 7, 5)
}

plot_interaction_model <- function(models, data, beyond, infield, file, model = "m02d") {
  m <- models[[which(vapply(models, `[[`, "", "name") == model)]]$fit
  q <- stats::quantile(data[[infield]], c(0.1, 0.5, 0.9), na.rm = TRUE)
  grid <- expand.grid(b = seq(min(data[[beyond]]), max(data[[beyond]]), length.out = 100), i = q)
  nd <- typical_rows(data, nrow(grid))
  nd[[beyond]] <- grid$b
  nd[[infield]] <- grid$i
  X <- stats::predict(m, nd, type = "lpmatrix")
  labs_s <- vapply(m$smooth, `[[`, "", "label")
  keep <- labs_s[grepl(beyond, labs_s, fixed = TRUE) | labs_s == paste0("s(", infield, ")")]
  cols <- unlist(lapply(m$smooth[labs_s %in% keep], function(s) s$first.para:s$last.para))
  eta <- drop(X[, cols] %*% stats::coef(m)[cols])
  se <- sqrt(rowSums((X[, cols] %*% m$Vp[cols, cols]) * X[, cols]))
  grid$level <- factor(grid$i, labels = c("low (10th percentile)", "median", "high (90th percentile)"))
  grid$fit <- exp(eta)
  grid$lwr <- exp(eta - 1.96 * se)
  grid$upr <- exp(eta + 1.96 * se)
  peaks <- dplyr::slice_max(grid, fit, n = 1, by = level)
  p <- ggplot(grid, aes(b, fit, colour = level, fill = level)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(data = peaks, size = 3, shape = 21, colour = "black") +
    geom_rug(data = data.frame(b = unique(data[[beyond]])), aes(b), inherit.aes = FALSE) +
    labs(x = paste("Beyond-field heterogeneity:", beyond), y = paste("Multiplicative effect on yield,", model),
         colour = "In-field complexity", fill = "In-field complexity") +
    theme_bw() +
    theme(legend.position = "bottom")
  save_plot(p, file, 7, 5)
}

significantterms <- function(models){
  df <- model_terms(models)
  df |> 
    dplyr::filter()
}