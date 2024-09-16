calulcate_distance <- function(ps, distance = "aitchison") {
  ps |>
    microViz::dist_calc(distance) |>
    update_provenance(ps, list(distance = distance))
}

test_distance <- function(ps, variable, limits) {
  loadNamespace("microViz")

  vi <- ps |> ps_variable_info(variable)

  if (!vi[["testable"]]) {
    cli::cli_alert_warning("{.var {variable}} must have at least two groups with at least two replicates each")
    return(invisible())
  }

  if (vi[[".length_levels"]] > limits[["variable_of_interest"]]) {
    cli::cli_alert_warning("{.var {variable}} must not have more than {.val {limits$variable_of_interest}} levels (has {.val {vi$.length_levels}})")
    return(invisible())
  }

  permanova <-
    ps |>
    microViz::dist_permanova(
      variables = variable,
      seed = 0L,
      verbose = FALSE
    ) |>
    microViz::perm_get()

  permanova_p_value <-
    permanova |>
    as.data.frame() |>
    tibble::rownames_to_column() |>
    dplyr::filter(rowname == variable) |>
    dplyr::pull("Pr(>F)")

  bdisp <-
    ps |>
    microViz::dist_bdisp(variables = variable) |>
    microViz::bdisp_get()

  bdisp_p_value <-
    bdisp |>
    purrr::chuck(variable) |>
    purrr::chuck("anova") |>
    as.data.frame() |>
    tibble::rownames_to_column() |>
    dplyr::filter(rowname == "Groups") |>
    dplyr::pull("Pr(>F)")

  ps |>
    provenance_as_tibble() |>
    tibble::add_column(
      `variable of interest` = variable,
      `PERMANOVA p-value` = permanova_p_value,
      `beta dispersion p-value` = bdisp_p_value
    )
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
