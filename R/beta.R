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

plot_ordination <- function(ps, variable, point_label, theme) {
  vi <- ps |> ps_variable_info(variable)

  if (!vi[["duplicates"]]) {
    cli::cli_alert_warning("{.var {variable}} must have duplicated values")
    return(invisible())
  }

  if (!vi[["multiple_levels"]]) {
    cli::cli_alert_warning("{.var {variable}} must have multiple levels")
    return(invisible())
  }

  plot <-
    ps |>
    microViz::ord_plot(
      # NOTE: color and colour behave differently in microViz
      colour = variable,
      fill = variable
    ) + ggplot2::labs(
      caption = NULL
    ) + theme + ggplot2::theme(
      aspect.ratio = 1.0
    )

  if (vi[["testable"]]) {
    plot <- plot + microViz::stat_chull(
      mapping = ggplot2::aes(colour = .data[[variable]], fill = .data[[variable]]),
      alpha = 0.1
    )
  }

  plot <- plot +
    ggrepel::geom_text_repel(
      mapping = ggplot2::aes(
        label = .data[[point_label]]
      ),
      max.overlaps = 100L,
      seed = 0L,
      size = 2.5,
      family = "Noto Sans"
    )

  plot |>
    plot_titles(
      ps,
      title = "beta diversity analysis"
    ) |>
    update_provenance(ps, list(
      aesthetics = list(color = variable)
    ))
}
