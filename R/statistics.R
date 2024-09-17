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
