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

plot_ordination <- function(ps, variable, point_label, limits, theme) {
  vi <- ps |> ps_variable_info(variable)
  vi_label <- ps |> ps_variable_info(point_label)

  if (!vi[["duplicates"]]) {
    cli::cli_alert_warning("{.var {variable}} must have duplicated values")
    return(invisible())
  }

  if (!vi[["multiple_levels"]]) {
    cli::cli_alert_warning("{.var {variable}} must have multiple levels")
    return(invisible())
  }

  if (vi[[".length_levels"]] > limits[["variable_of_interest"]]) {
    cli::cli_alert_warning("{.var {variable}} must not have more than {.val {limits$variable_of_interest}} levels (has {.val {vi$.length_levels}})")
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

  if (vi_label[["multiple_levels"]] && vi[[".length"]] <= limits[["sample"]]) {
    plot <- plot +
      ggrepel::geom_text_repel(
        mapping = ggplot2::aes(
          label = .data[[point_label]]
        ),
        max.overlaps = 100L,
        seed = 0L,
        size = 2.5,
        family = font_family()
      )
  }

  plot |>
    plot_titles(
      ps,
      title = "beta diversity analysis"
    ) |>
    update_provenance(ps, list(
      aesthetics = list(color = variable)
    ))
}
