calulcate_distance <- function(ps, distance = "aitchison") {
  if (is.null(ps)) {
    return(invisible())
  }

  ps |>
    microViz::dist_calc(distance) |>
    update_provenance(ps, list(distance = distance))
}

calulcate_ordination <- function(ps) {
  if (is.null(ps)) {
    return(invisible())
  }

  loadNamespace("microViz")

  if (phyloseq::nsamples(ps) < 2L) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(ps)}}: at least two samples are required for ordination analysis")
    return(invisible())
  }

  ps_new <-
    ps |>
    microViz::ord_calc()

  ps_new |>
    update_provenance(ps, list(
      ordination = microViz::info_get(ps_new)[["ord_info"]][["method"]]
    ))
}

format_beta_diversity_test <- function(beta_diversity_test_raw) {
  if (is.null(beta_diversity_test_raw)) {
    return(invisible())
  }

  beta_diversity_test_raw |>
    dplyr::rename(metric = distance)
}

plot_ordination <- function(ps, variable = "Group", point_label, limits, theme = ggplot_theme()) {
  if (is.null(ps)) {
    return(invisible())
  }

  # this is in line with test_distance()
  ps <- ps |> microViz::ps_mutate(across(any_of(variable), as.factor))

  vi <- ps |> ps_variable_info(variable)
  vi_label <- ps |> ps_variable_info(point_label)

  if (!vi[["duplicates"]]) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(ps)}}: {.var {variable}} must have duplicated values")
    return(invisible())
  }

  if (!vi[["multiple_levels"]]) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(ps)}}: {.var {variable}} must have multiple levels")
    return(invisible())
  }

  if (vi[[".length_levels"]] > limits[["variable_of_interest"]]) {
    cli::cli_alert_warning(str_c(
      "{.field {provenance_as_short_title(ps)}}: ",
      "{.var {variable}} must not have more than {.val {limits$variable_of_interest}} levels (has {.val {vi$.length_levels}})"
    ))
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
    update_provenance(ps, list(
      aesthetics = list(color = variable)
    ))
}

plot_ordination_with_tests <- function(plot, test_result) {
  if (is.null(plot)) {
    return(invisible())
  }

  variable_of_interest <-
    test_result |>
    dplyr::pull(`variable of interest`) |>
    unique()

  # assertion that plot and test belong together
  stopifnot(identical(
    plot |> get_provenance() |> purrr::list_assign(ordination = rlang::zap()),
    test_result |> get_provenance() |> purrr::list_assign(aesthetics = list(color = variable_of_interest))
  ))

  permanova_p_value <-
    test_result |>
    dplyr::filter(test == "PERMANOVA", NAs != "removed") |>
    dplyr::pull(`p-value`)
  beta_dispersion_p_value <-
    test_result |>
    dplyr::filter(test == "beta dispersion ANOVA", NAs != "removed") |>
    dplyr::pull(`p-value`)

  subtitle <- ""
  if (!is.na(permanova_p_value)) {
    p_format(0L) # does nothing, just to enable re-evaluation
    subtitle <- glue::glue("{variable_of_interest} test: PERMANOVA {p_format(permanova_p_value)}")
    if (!is.na(beta_dispersion_p_value)) {
      subtitle <- subtitle |> stringr::str_c(glue::glue(" with dispersion ANOVA {p_format(beta_dispersion_p_value)}"))
    }
  }

  permanova_p_value <-
    test_result |>
    dplyr::filter(test == "PERMANOVA", NAs == "removed") |>
    dplyr::pull(`p-value`)
  beta_dispersion_p_value <-
    test_result |>
    dplyr::filter(test == "beta dispersion ANOVA", NAs == "removed") |>
    dplyr::pull(`p-value`)

  if (length(permanova_p_value) > 0L && !is.na(permanova_p_value)) {
    subtitle <- subtitle |> stringr::str_c("\n", glue::glue("\u21b3without NAs: PERMANOVA {p_format(permanova_p_value)}"))
    if (length(beta_dispersion_p_value) > 0L && !is.na(beta_dispersion_p_value)) {
      subtitle <- subtitle |> stringr::str_c(glue::glue(" with dispersion ANOVA {p_format(beta_dispersion_p_value)}"))
    }
  }

  title <- str_c(
    "beta diversity analysis (",
    plot |> get_provenance() |> chuck("distance"), " ",
    plot |> get_provenance() |> chuck("ordination"), ")"
  )

  plot |>
    plot_titles(
      title = title,
      subtitles = subtitle,
      distance = zap(),
      ordination = zap()
    )
}
