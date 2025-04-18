p_format <- function(x) {
  withr::local_options(list(scipen = 999L))

  x |>
    rstatix::p_format(add.p = TRUE) |>
    rstatix::p_mark_significant() |>
    dplyr::na_if("p=NA")
}

finalize_tests_table <- function(data) {
  data |>
    dplyr::mutate(`p-value formatted` = `p-value` |> p_format(), .after = `p-value`)
}

super_safely <- function(fun, ...) {
  res <- purrr::quietly(purrr::safely(fun))(...) |> purrr::list_flatten(name_spec = "{inner}")

  res[["log"]] <-
    res[rev(c("messages", "warnings", "error"))] |>
    purrr::map(as.character) |>
    unlist() |>
    stringr::str_trim() |>
    stringr::str_flatten(" | ")

  res
}

test_distance <- function(ps, variable) {
  if (is.null(ps)) {
    cli::cli_alert_warning("skipped because container is NULL")
    return(invisible())
  }

  loadNamespace("microViz")

  # all variables are tested as categories (factors) - even continuous variables (for now)
  ps <- ps |> microViz::ps_mutate(across(any_of(variable), fortify))

  permanova_p_value <- NA_real_
  permanova_error <- NULL
  bdisp_p_value <- NA_real_
  bdisp_error <- NULL

  vi <- ps |> ps_variable_info(variable)

  if (!vi[["testable"]]) { # nolint: if_not_else_linter.
    permanova_error <- bdisp_error <- "must have at least two groups with at least two replicates each"
  } else {
    permanova <- super_safely(
      microViz::dist_permanova,
      ps,
      variables = variable,
      n_perms = 9999L,
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

  seq_len(2L) |>
    map(\(i) provenance_as_tibble(ps)) |>
    dplyr::bind_rows() |>
    tibble::add_column(
      `variable of interest` = variable,
      test = c("PERMANOVA", "beta dispersion ANOVA"),
      `p-value` = c(permanova_p_value, bdisp_p_value),
      error = c(permanova_error |> cli::pluralize(), bdisp_error |> cli::pluralize())
    ) |>
    update_provenance(ps)
}
