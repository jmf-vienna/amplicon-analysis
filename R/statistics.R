test_distance <- function(ps, variable, limits) {
  loadNamespace("microViz")

  permanova_p_value <- NA_real_
  permanova_error <- NULL
  bdisp_p_value <- NA_real_
  bdisp_error <- NULL

  vi <- ps |> ps_variable_info(variable)

  if (!vi[["testable"]]) {
    permanova_error <- "must have at least two groups with at least two replicates each"
  } else if (vi[[".length_levels"]] > limits[["variable_of_interest"]]) {
    permanova_error <- "must not have more than {limits$variable_of_interest} levels (has {vi$.length_levels})"
  }

  if (is.null(permanova_error)) {
    ps_extra <-
      ps |>
      microViz::dist_permanova(
        variables = variable,
        seed = 0L,
        verbose = FALSE
      )

    permanova_p_value <-
      ps_extra |>
      microViz::perm_get() |>
      as.data.frame() |>
      tibble::rownames_to_column() |>
      dplyr::filter(rowname == variable) |>
      dplyr::pull("Pr(>F)")

    ps_extra <-
      ps_extra |>
      microViz::dist_bdisp(variables = variable)

    bdisp_p_value <-
      ps_extra |>
      microViz::bdisp_get() |>
      purrr::chuck(variable) |>
      purrr::chuck("anova") |>
      as.data.frame() |>
      tibble::rownames_to_column() |>
      dplyr::filter(rowname == "Groups") |>
      dplyr::pull("Pr(>F)")
  }

  if (is.null(permanova_error)) {
    permanova_error <- ""
  } else {
    cli::cli_alert_warning(stringr::str_c("{.var {variable}} ", permanova_error))
  }

  if (is.null(bdisp_error)) {
    bdisp_error <- ""
  } else {
    cli::cli_alert_warning(stringr::str_c("{.var {variable}} ", bdisp_error))
  }

  ps |>
    provenance_as_tibble() |>
    tibble::add_column(
      `variable of interest` = variable,
      `PERMANOVA p-value` = permanova_p_value,
      `PERMANOVA error` = permanova_error |> cli::pluralize(),
      `beta dispersion p-value` = bdisp_p_value,
      `beta dispersion error` = bdisp_error |> cli::pluralize()
    ) |>
    update_provenance(ps)
}
