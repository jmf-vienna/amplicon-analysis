p_format <- function(x) {
  withr::local_options(list(scipen = 999L))

  x |>
    rstatix::p_format(add.p = TRUE) |>
    rstatix::p_mark_significant() |>
    dplyr::na_if("p=NA")
}

single_test <- function(x, two_sample_test, variable, value) {
  x <- dplyr::mutate(x, across(any_of(variable), droplevels))

  test_formula <- as.formula(str_c(value, "~", variable))

  if (identical(two_sample_test, "wilcox")) {
    test_result <-
      x |>
      wilcox.test(test_formula, data = _) |>
      # "cannot compute exact p-value with ties" warnings are ignored
      suppressWarnings()
  }

  tibble(
    .y. = value,
    group1 = x |> pull(variable) |> levels() |> dplyr::first(),
    group2 = x |> pull(variable) |> levels() |> dplyr::last(),
    p = test_result |> chuck("p.value")
  )
}

pairwise_test <- function(data, variable, value, two_sample_test, p_adjust_method) {
  data |>
    purrr::chuck(variable) |>
    jmf::uniques() |>
    utils::combn(2L, simplify = FALSE) |>
    purrr::map(\(x) {
      data |>
        dplyr::filter(.data[[variable]] %in% x) |>
        single_test(two_sample_test, variable, value)
    }) |>
    dplyr::bind_rows() |>
    rstatix::adjust_pvalue(method = p_adjust_method) |>
    rstatix::add_significance("p.adj")
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

test_distance <- function(ps_raw, variable, .filter_na = FALSE) { # nolint: cyclocomp_linter.
  if (is.null(ps_raw)) {
    return(invisible())
  }

  vi_raw <- ps_variable_info(ps_raw, variable)

  if (!vi_raw[["exists"]]) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(ps_raw)}}: skipped because variable does not exist")
    return(invisible())
  }

  loadNamespace("microViz")

  if (.filter_na) {
    # this filters only the phyloseq object, but not the distance!
    # dist_permanova and bdisp_get then call ps_drop_incomplete (again) and filter @dist
    # so this works, but could have unexpected side effects
    ps_raw <- ps_raw |> microViz::ps_drop_incomplete(variable)
  }

  # all variables are tested as categories (factors) - even continuous variables (for now)

  ps <- ps_raw |> microViz::ps_mutate(across(any_of(variable), fortify))
  vi <- ps_variable_info(ps, variable)

  permanova_p_value <- NA_real_
  permanova_error <- NULL
  bdisp_p_value <- NA_real_
  bdisp_error <- NULL

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
    cli::cli_alert_warning(stringr::str_c("{.field {provenance_as_short_title(ps_raw)}}: {.var {variable}} PERMANOVA: ", permanova_error))
  }

  if (is.null(bdisp_error)) {
    bdisp_error <- ""
  } else {
    cli::cli_alert_warning(stringr::str_c("{.field {provenance_as_short_title(ps_raw)}}: {.var {variable}} beta dispertion: ", bdisp_error))
  }

  res <-
    seq_len(2L) |>
    map(\(i) provenance_as_tibble(ps)) |>
    dplyr::bind_rows() |>
    tibble::add_column(
      `variable of interest` = variable,
      NAs = case_when(
        !vi_raw[["has_na"]] ~ "none",
        vi_raw[["all_na"]] ~ "all",
        .filter_na ~ "removed",
        !.filter_na ~ "enumerated"
      ),
      test = c("PERMANOVA", "beta dispersion ANOVA"),
      `p-value` = c(permanova_p_value, bdisp_p_value),
      error = c(permanova_error |> cli::pluralize(), bdisp_error |> cli::pluralize())
    )

  if (!.filter_na && vi_raw[["has_na"]] && !vi_raw[["all_na"]]) {
    res <- res |> bind_rows(test_distance(ps_raw, variable, TRUE))
  }

  res |> update_provenance(ps)
}
