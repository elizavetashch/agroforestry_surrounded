# ---- helpers ---------------------------------------------------------------

control_field_colours <- col_system

summarise_yield <- function(d, groups) {
  d |>
    dplyr::filter(!is.na(yield_tha)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
    dplyr::summarise(
      mean = mean(yield_tha),
      n    = dplyr::n(),
      se   = dplyr::if_else(n > 1, stats::sd(yield_tha) / sqrt(n), NA_real_),
      .groups = "drop"
    )
}

# keep only field-year-crop combinations that have both control and agroforestry
keep_paired <- function(s) {
  s |>
    dplyr::group_by(field, year, crop_unified) |>
    dplyr::filter(dplyr::n_distinct(control_field) == 2) |>
    dplyr::ungroup()
}



# ---- (3) all fields in one graph: control yield ---------------------------

control_yield <- function(d) {
  summarise_yield(d, c("field", "year", "crop_unified", "control_field")) |>
    keep_paired() |>
    dplyr::select(field, year, crop_unified, control_field, mean) |>
    tidyr::pivot_wider(names_from = control_field, values_from = mean) |>
    dplyr::mutate(rel_yield = agroforestry / control * 100)
}

plot_control_yield <- function(rel) {
  ggplot2::ggplot(rel, ggplot2::aes(field, rel_yield)) +
    ggplot2::geom_hline(yintercept = 100, linetype = "dashed") +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey90") +
    ggplot2::geom_jitter(ggplot2::aes(colour = crop_unified),
                         width = 0.15, height = 0, size = 2) +
    ggsci::scale_colour_npg() +
    ggplot2::labs(x = NULL,
                  y = "Agroforestry yield control to control (%)",
                  colour = "Crop") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


plot_yield_timeseries <- function(d, file) {
  s <- summarise_yield(d, c("field", "year", "crop_unified", "control_field"))
  
  crops <- s |>
    dplyr::group_by(field, year) |>
    dplyr::summarise(crop = paste(sort(unique(crop_unified)), collapse = "/"),
                     .groups = "drop")
  
  pd <- position_dodge(width = 0.3)
  
  p <- ggplot(s, aes(year, mean, colour = control_field, group = control_field)) +
    geom_line(position = pd) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.2, position = pd, na.rm = TRUE) +
    geom_point(size = 2, position = pd) +
    geom_text(data = crops, aes(x = year, y = -Inf, label = crop),
              inherit.aes = FALSE, angle = 90, hjust = -0.1, size = 3) +
    facet_wrap(~ field, scales = "free") +
    scale_colour_manual(values = control_field_colours) +
    scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)) +
    labs(x = "Year", y = "Yield (t/ha), mean ± SE", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14)
    )
  print(p)
  save_plot(p, file, 14, 9)
}

plot_yield_bars <- function(d, file) {
  s <- summarise_yield(d, c("field", "year", "crop_unified", "control_field")) |>
    keep_paired() |>
    dplyr::mutate(label = paste0(year, "\n", crop_unified))
  
  dodge <- position_dodge(width = 0.8)
  
  p <- ggplot(s, aes(label, mean, fill = control_field)) +
    geom_col(position = dodge, width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = dodge, width = 0.2, na.rm = TRUE) +
    facet_wrap(~ field, scales = "free") +
    scale_fill_manual(values = control_field_colours) +
    labs(x = NULL, y = "Yield (t/ha), mean ± SE", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14)
    )
  print(p)
  save_plot(p, file, 14, 9)
}

plot_control_yield <- function(rel, file) {
  p <- ggplot(rel, aes(field, rel_yield)) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
    geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_jitter(aes(colour = crop_unified), width = 0.15, height = 0, size = 2) +
    ggsci::scale_colour_npg() +
    labs(x = NULL, y = "Agroforestry yield control to control (%)", colour = "Crop") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  print(p)
  save_plot(p, file, 10, 7)
}