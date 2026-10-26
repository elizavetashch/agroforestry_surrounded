# Data exploration protocol after Zuur, Ieno & Elphick (2010),
# Methods Ecol Evol 1: 3-14. One function per step; each saves its figures and
# tables to `dir` and returns the file paths.

numeric_cols <- function(df, cols) cols[vapply(df[cols], is.numeric, logical(1))]

# Step 1: outliers in Y and X (Cleveland dotplots)
zuur_step1_outliers <- function(df, y, covariates, dir) {
  vars <- numeric_cols(df, c(y, covariates))
  long <- df |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::select(.row, dplyr::all_of(vars)) |>
    tidyr::pivot_longer(-.row, names_to = "variable", values_to = "value")

  p <- ggplot(long, aes(value, .row)) +
    geom_point(colour = "steelblue", size = 0.8) +
    facet_wrap(~ variable, scales = "free_x") +
    labs(x = "Value", y = "Order of the data") +
    theme_bw()

  tab <- long |>
    dplyr::summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      n_beyond_3sd = sum(abs(value - mean) > 3 * sd, na.rm = TRUE),
      .by = variable
    )

  c(save_plot(p, file.path(dir, "step1_dotplots.png"), 10, 8),
    save_table(tab, file.path(dir, "step1_outlier_summary.csv")))
}

# Step 2: homogeneity of variance (conditional boxplots, Levene test centred on the median)
levene <- function(y, g) {
  g <- factor(g)
  dev <- abs(y - stats::ave(y, g, FUN = function(x) stats::median(x, na.rm = TRUE)))
  a <- stats::anova(stats::lm(dev ~ g))
  c(df1 = a$Df[1], df2 = a$Df[2], F = a$`F value`[1], p_value = a$`Pr(>F)`[1])
}

zuur_step2_variance <- function(df, y, groups, dir) {
  plots <- vapply(groups, function(g) {
    p <- ggplot(df, aes(factor(.data[[g]]), .data[[y]])) +
      geom_boxplot(fill = "grey85", outlier.colour = "red") +
      labs(x = g, y = y) +
      theme_bw()
    save_plot(p, file.path(dir, paste0("step2_boxplot_", g, ".png")))
  }, character(1))

  tab <- dplyr::bind_rows(lapply(groups, function(g) {
    data.frame(group = g, t(levene(df[[y]], df[[g]])))
  }))

  c(plots, save_table(tab, file.path(dir, "step2_levene.csv")))
}

# Step 3: normality of Y, overall and within groups
shapiro <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(x) > 5000 || length(unique(x)) < 2) return(c(W = NA, p_value = NA))
  t <- stats::shapiro.test(x)
  c(W = unname(t$statistic), p_value = t$p.value)
}

zuur_step3_normality <- function(df, y, groups, dir) {
  p_hist <- ggplot(df, aes(.data[[y]])) +
    geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
    theme_bw()
  p_qq <- ggplot(df, aes(sample = .data[[y]])) +
    stat_qq(colour = "steelblue") +
    stat_qq_line(colour = "red") +
    labs(x = "Theoretical quantiles", y = y) +
    theme_bw()

  plots <- vapply(groups, function(g) {
    p <- ggplot(df, aes(.data[[y]])) +
      geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
      facet_wrap(vars(.data[[g]]), scales = "free_y") +
      theme_bw()
    save_plot(p, file.path(dir, paste0("step3_histogram_by_", g, ".png")), 10, 7)
  }, character(1))

  tab <- dplyr::bind_rows(
    data.frame(group = "all", level = "all", n = sum(!is.na(df[[y]])), t(shapiro(df[[y]]))),
    lapply(groups, function(g) {
      dplyr::bind_rows(lapply(split(df[[y]], df[[g]]), function(x) {
        data.frame(n = sum(!is.na(x)), t(shapiro(x)))
      }), .id = "level") |>
        dplyr::mutate(group = g, .before = 1)
    })
  )

  c(save_plot(p_hist, file.path(dir, "step3_histogram.png")),
    save_plot(p_qq, file.path(dir, "step3_qqplot.png")),
    plots,
    save_table(tab, file.path(dir, "step3_shapiro.csv")))
}

# Step 4: zeros in Y
zuur_step4_zeros <- function(df, y, groups, dir) {
  tab <- dplyr::bind_rows(
    data.frame(group = "all", level = "all", n = nrow(df), n_zero = sum(df[[y]] == 0, na.rm = TRUE)),
    lapply(groups, function(g) {
      df |>
        dplyr::summarise(n = dplyr::n(), n_zero = sum(.data[[y]] == 0, na.rm = TRUE), .by = dplyr::all_of(g)) |>
        dplyr::transmute(group = g, level = as.character(.data[[g]]), n, n_zero)
    })
  ) |>
    dplyr::mutate(pct_zero = round(100 * n_zero / n, 1))

  save_table(tab, file.path(dir, "step4_zeros.csv"))
}

# Step 5: collinearity among covariates
vif <- function(df, vars) {
  d <- stats::na.omit(df[vars])
  vapply(vars, function(v) {
    r2 <- summary(stats::lm(stats::reformulate(setdiff(vars, v), v), data = d))$r.squared
    1 / (1 - r2)
  }, numeric(1))
}

zuur_step5_collinearity <- function(df, covariates, dir) {
  vars <- numeric_cols(df, covariates)
  cor_tab <- round(stats::cor(df[vars], use = "pairwise.complete.obs"), 2)
  cor_tab <- data.frame(variable = rownames(cor_tab), cor_tab, row.names = NULL)
  vif_tab <- data.frame(variable = vars, vif = round(vif(df, vars), 2), row.names = NULL)

  f_pairs <- file.path(dir, "step5_pairs.png")
  grDevices::png(f_pairs, width = 3000, height = 3000, res = 200)
  graphics::pairs(df[vars], panel = graphics::panel.smooth, pch = 16, cex = 0.4, col = "grey40")
  grDevices::dev.off()

  c(f_pairs,
    save_table(cor_tab, file.path(dir, "step5_correlation.csv")),
    save_table(vif_tab, file.path(dir, "step5_vif.csv")))
}

# Step 6: relationships Y ~ X, conditioned on a grouping factor
zuur_step6_relationships <- function(df, y, covariates, condition, dir) {
  vapply(numeric_cols(df, covariates), function(x) {
    p <- ggplot(df, aes(.data[[x]], .data[[y]])) +
      geom_point(alpha = 0.4, size = 1) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "steelblue") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "red", linetype = "dashed") +
      facet_wrap(vars(.data[[condition]])) +
      labs(x = x, y = y) +
      theme_bw()
    save_plot(p, file.path(dir, paste0("step6_", y, "_vs_", x, ".png")))
  }, character(1))
}

# Step 7: interactions (non-parallel lines suggest an interaction)
zuur_step7_interactions <- function(df, y, covariates, groups, dir) {
  combos <- expand.grid(x = numeric_cols(df, covariates), g = groups, stringsAsFactors = FALSE)
  mapply(function(x, g) {
    p <- ggplot(df, aes(.data[[x]], .data[[y]], colour = factor(.data[[g]]), fill = factor(.data[[g]]))) +
      geom_smooth(method = "lm", formula = y ~ x, alpha = 0.15) +
      labs(x = x, y = y, colour = g, fill = g) +
      theme_bw()
    save_plot(p, file.path(dir, paste0("step7_", x, "_x_", g, ".png")))
  }, combos$x, combos$g, USE.NAMES = FALSE)
}

# Step 8: independence (residuals of a naive model by group, ACF, variance components)
zuur_step8_independence <- function(df, naive_formula, random, order_by, dir) {
  m <- stats::lm(stats::as.formula(naive_formula), data = df, na.action = stats::na.exclude)
  d <- dplyr::mutate(df, .resid = stats::resid(m))

  plots <- vapply(random, function(g) {
    p <- ggplot(d, aes(factor(.data[[g]]), .resid)) +
      geom_boxplot(fill = "grey85") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      labs(x = g, y = "Residual (naive model)") +
      theme_bw()
    save_plot(p, file.path(dir, paste0("step8_residuals_by_", g, ".png")))
  }, character(1))

  res <- stats::na.omit(dplyr::arrange(d, dplyr::across(dplyr::all_of(order_by)))$.resid)
  a <- stats::acf(res, plot = FALSE)
  ci <- stats::qnorm(0.975) / sqrt(length(res))
  p_acf <- ggplot(data.frame(lag = drop(a$lag), acf = drop(a$acf)), aes(lag, acf)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", colour = "blue") +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(x = "Lag", y = "ACF of naive residuals") +
    theme_bw()

  re_terms <- paste0("s(", random, ", bs = \"re\")", collapse = " + ")
  dd <- dplyr::mutate(df, dplyr::across(dplyr::all_of(random), as.factor))
  g <- mgcv::gam(stats::as.formula(paste(naive_formula, "+", re_terms)), data = dd, method = "REML")
  utils::capture.output(vc <- mgcv::gam.vcomp(g, rescale = FALSE))
  
  sd <- if (is.matrix(vc)) vc[, "std.dev"] else if (is.list(vc)) vc$all else vc
  labs <- c(vapply(g$smooth, function(s) s$label, ""), "scale")
  comp <- if (!is.null(names(sd))) names(sd) else if (length(labs) == length(sd)) labs else paste0("component_", seq_along(sd))
  vc_tab <- data.frame(component = comp, variance = unname(sd)^2, row.names = NULL) |>
    dplyr::mutate(pct_variance = round(100 * variance / sum(variance), 1))
  
  c(plots,
    save_plot(p_acf, file.path(dir, "step8_acf.png")),
    save_table(vc_tab, file.path(dir, "step8_variance_components.csv")))
}



# Further Investigations behind Zuur --------------------------------------


# (1) What crops are winter crops what crops are summer crops  ------------

season_croptype <- function(df, dir){
  seasontab <- table(df$crop_unified, df$crop_season)
  save_table(seasontab, file.path(dir, "add_season_croptype.csv"))
}


# (2) Distribution choice -------------------------------------------------

# Step 9: which distribution fits the response?

family_qq <- function(df, y, dir) {
  x <- df[[y]]
  x <- x[!is.na(x) & x > 0]
  p <- stats::ppoints(length(x))
  m <- mean(x)
  g <- MASS::fitdistr(x / m, "gamma")$estimate
  ln <- MASS::fitdistr(x, "lognormal")$estimate
  w <- suppressWarnings(MASS::fitdistr(x / m, "weibull")$estimate)
  dists <- list(
    Normal    = list(q = stats::qnorm(p, m, stats::sd(x)), ll = sum(stats::dnorm(x, m, stats::sd(x), log = TRUE))),
    Lognormal = list(q = stats::qlnorm(p, ln[1], ln[2]), ll = sum(stats::dlnorm(x, ln[1], ln[2], log = TRUE))),
    Gamma     = list(q = stats::qgamma(p, g[1], g[2] / m), ll = sum(stats::dgamma(x, g[1], g[2] / m, log = TRUE))),
    Weibull   = list(q = stats::qweibull(p, w[1], w[2] * m), ll = sum(stats::dweibull(x, w[1], w[2] * m, log = TRUE)))
  )
  qq <- dplyr::bind_rows(lapply(names(dists), function(d) {
    data.frame(distribution = d, theoretical = dists[[d]]$q, observed = sort(x))
  }))
  p_qq <- ggplot(qq, aes(theoretical, observed)) +
    geom_abline(colour = "red") +
    geom_point(size = 0.8, alpha = 0.5) +
    facet_wrap(~ distribution, scales = "free") +
    labs(x = "Theoretical quantiles", y = paste("Observed", y)) +
    theme_bw()
  tab <- data.frame(distribution = names(dists), logLik = vapply(dists, `[[`, 0, "ll")) |>
    dplyr::mutate(AIC = -2 * logLik + 2 * 2, delta_AIC = AIC - min(AIC)) |>
    dplyr::arrange(AIC)
  c(save_plot(p_qq, file.path(dir, "step9_qq_distributions.png"), 8, 7),
    save_table(tab, file.path(dir, "step9_distribution_fit.csv")))
}

family_compare <- function(data, formula, dir) {
  f <- stats::as.formula(formula)
  y <- data[[all.vars(f)[1]]]
  fits <- list(
    Gaussian  = mgcv::gam(f, data = data, family = gaussian(), method = "ML"),
    Gamma_log = mgcv::gam(f, data = data, family = Gamma(link = "log"), method = "ML"),
    Tweedie   = mgcv::gam(f, data = data, family = mgcv::tw(), method = "ML"),
    Lognormal = mgcv::gam(stats::update(f, log(.) ~ .), data = data, family = gaussian(), method = "ML")
  )
  tab <- dplyr::bind_rows(lapply(names(fits), function(n) {
    a <- stats::AIC(fits[[n]])
    if (n == "Lognormal") a <- a + 2 * sum(log(y))   # puts the log-scale model on the yield scale
    data.frame(family = n, AIC = a, deviance_explained = summary(fits[[n]])$dev.expl)
  })) |>
    dplyr::mutate(delta_AIC = AIC - min(AIC), .after = AIC) |>
    dplyr::arrange(AIC)
  
  f_png <- file.path(dir, "family_residuals.png")
  grDevices::png(f_png, width = 4000, height = 2000, res = 300)
  graphics::par(mfrow = c(2, 4))
  for (n in names(fits)) mgcv::qq.gam(fits[[n]], pch = 16, cex = 0.4, main = paste("QQ:", n))
  for (n in names(fits)) {
    plot(stats::fitted(fits[[n]]), stats::residuals(fits[[n]], type = "deviance"), pch = 16, cex = 0.4,
         col = "grey40", xlab = "Fitted values", ylab = "Deviance residuals", main = paste("Residuals:", n))
    graphics::abline(h = 0, col = "red", lty = 2)
  }
  grDevices::dev.off()
  
  c(save_table(tab, file.path(dir, "step9_family_comparison_aic.csv")), f_png)
}


mean_sd_compare <- function(df, dir){
    meansdtab <- 
      df |> 
      dplyr::summarise(n = dplyr::n(), 
                     mean = mean(yield_tha), 
                     sd = sd(yield_tha),
                     cv = sd / mean, .by = crop_unified) |>
    dplyr::arrange(mean)
  save_table(meansdtab, file.path(dir, "step9_mean_sd_compare.csv"))
}
