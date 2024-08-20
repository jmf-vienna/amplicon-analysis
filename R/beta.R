calulcate_distance <- function(ps) {
  ps |>
    microViz::dist_calc("aitchison")
}

calulcate_ordination <- function(ps) {
  ps |>
    microViz::ord_calc()
}

plot_ordination <- function(ps, group) {
  ps |>
    microViz::ord_plot(
      colour = group,
      fill = group
    ) + ggplot2::theme_get()
}
