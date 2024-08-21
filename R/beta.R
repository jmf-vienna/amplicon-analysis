calulcate_distance <- function(ps) {
  ps |>
    microViz::dist_calc("aitchison") |>
    set_provenance(ps)
}

calulcate_ordination <- function(ps) {
  ps |>
    microViz::ord_calc() |>
    set_provenance(ps)
}

plot_ordination <- function(ps, group, theme) {
  plot <-
    ps |>
    microViz::ord_plot(
      colour = group,
      fill = group
    ) + ggplot2::labs(
      title = ps |> collapse_provenance(" | ")
    ) + theme
  plot |> set_provenance(ps, list(plot = "ordination"))
}
