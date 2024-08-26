calulcate_distance <- function(ps, distance = "aitchison") {
  ps |>
    microViz::dist_calc(distance) |>
    update_provenance(ps, list(distance = distance))
}

calulcate_ordination <- function(ps) {
  ps_new <-
    ps |>
    microViz::ord_calc()

  ps_new |>
    update_provenance(ps, list(
      ordination = microViz::info_get(ps_new)[["ord_info"]][["method"]]
    ))
}

plot_ordination <- function(ps, group, theme) {
  plot <-
    ps |>
    microViz::ord_plot(
      color = group,
      fill = group
    ) + ggplot2::labs(
      title = ps |> as_title()
    ) + theme

  plot |> update_provenance(ps, list(
    aesthetics = list(color_by = group)
  ))
}
