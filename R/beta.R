calulcate_distance <- function(ps) {
  ps |>
    microViz::dist_calc("aitchison") |>
    update_provenance(ps)
}

calulcate_ordination <- function(ps) {
  ps |>
    microViz::ord_calc() |>
    update_provenance(ps)
}

plot_ordination <- function(ps, group, theme) {
  plot <-
    ps |>
    microViz::ord_plot(
      colour = group,
      fill = group
    ) + ggplot2::labs(
      title = ps |> as_title()
    ) + theme

  plot |> update_provenance(ps, list(
    plot = "ordination",
    aesthetics = list(color_by = group)
  ))
}
