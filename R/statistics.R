super_safely <- function(fun, ...) {
  res <- purrr::quietly(purrr::safely(fun))(...) |> purrr::list_flatten(name_spec = "{inner}")

  res[["log"]] <-
    res[c("messages", "warnings", "error")] |>
    purrr::map(as.character) |>
    unlist() |>
    stringr::str_trim() |>
    stringr::str_flatten(" | ")

  res
}

test_distance <- function(ps, variable, limits) {
  loadNamespace("microViz")

  permanova_p_value <- NA_real_
  permanova_error <- NULL
  bdisp_p_value <- NA_real_
  bdisp_error <- NULL

  vi <- ps |> ps_variable_info(variable)

  if (!vi[["testable"]]) {
    permanova_error <- bdisp_error <- "must have at least two groups with at least two replicates each"
  } else {
    permanova <- super_safely(
      microViz::dist_permanova,
      ps,
      variables = variable,
      seed = 0L,
      verbose = FALSE
    )

    if (permanova[["log"]] == "") {
      permanova_p_value <-
        permanova |>
        purrr::chuck("result") |>
        microViz::perm_get() |>
        as.data.frame() |>
        tibble::rownames_to_column() |>
        dplyr::filter(rowname == variable) |>
        dplyr::pull("Pr(>F)")
    } else {
      permanova_error <- permanova[["log"]]
    }

    bdisp <- super_safely(
      microViz::dist_bdisp,
      ps,
      variables = variable
    )

    if (bdisp[["log"]] == "") {
      bdisp_p_value <-
        bdisp |>
        purrr::chuck("result") |>
        microViz::bdisp_get() |>
        purrr::chuck(variable) |>
        purrr::chuck("anova") |>
        as.data.frame() |>
        tibble::rownames_to_column() |>
        dplyr::filter(rowname == "Groups") |>
        dplyr::pull("Pr(>F)")
    } else {
      bdisp_error <- bdisp[["log"]]
    }
  }

  if (is.null(permanova_error)) {
    permanova_error <- ""
  } else {
    cli::cli_alert_warning(stringr::str_c("{.var {variable}} PERMANOVA: ", permanova_error))
  }

  if (is.null(bdisp_error)) {
    bdisp_error <- ""
  } else {
    cli::cli_alert_warning(stringr::str_c("{.var {variable}} beta dispertion: ", bdisp_error))
  }

  ps |>
    provenance_as_tibble() |>
    tibble::add_column(
      `variable of interest` = variable,
      `PERMANOVA p-value` = permanova_p_value,
      `PERMANOVA error` = permanova_error |> cli::pluralize(),
      `beta dispersion ANOVA p-value` = bdisp_p_value,
      `beta dispersion error` = bdisp_error |> cli::pluralize()
    ) |>
    update_provenance(ps)
}
