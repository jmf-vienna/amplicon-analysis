add_alpha_diversity <- function(se, alpha_diversity_indexes = "observed", threshold = 1000L, rarefaction_rounds = 10L) {
  min <- min_col_sum(se)
  if (min >= threshold) {
    se <-
      se |>
      mia::addAlpha(
        index = alpha_diversity_indexes,
        name = str_c(".alpha_diversity_", alpha_diversity_indexes)
      )

    # rarefaction makes no sense if there is only one sample
    if (ncol(se) > 1L) {
      set.seed(0L)
      se <-
        se |>
        mia::addAlpha(
          index = alpha_diversity_indexes,
          name = str_c(".alpha_diversity_at_", min, "_", alpha_diversity_indexes),
          sample = min,
          niter = rarefaction_rounds
        ) |>
        # suppress vegan::rrarefy warning message "function should be used for observed counts, but smallest count is %d"
        suppressWarnings()
    }
  } else {
    cli::cli_alert_warning(
      "{.field {provenance_as_short_title(se)}}: alpha diversity skipped (minimum sample size {.val {min}} is below threshold {.val {threshold}})"
    )
  }
  se
}

get_alpha_diversity <- function(se, variable, theme) {
  res <-
    se |>
    SummarizedExperiment::colData() |>
    as_full_tibble("ID")

  if (res |> select(starts_with(".alpha_diversity_")) |> ncol() |> identical(0L)) {
    # return empty tibble if no alpha diversity indexes are present
    res <-
      res |>
      add_column(Index = NA, Rarefaction = NA, Diversity = NA) |>
      dplyr::filter(is.na(ID))
  } else {
    res <-
      res |>
      pivot_longer(starts_with(".alpha_diversity_"), names_to = "Index", values_to = "Diversity") |>
      dplyr::filter(!str_ends(Index, fixed("_se"))) |>
      mutate(
        Rarefaction = Index |> str_extract("at_[0-9]+") |> str_replace(fixed("at_"), "rarefaction depth ") |> str_replace_na("no rarefaction"),
        Index = Index |> str_remove("^[.]alpha_diversity_(at_[0-9]+_)?") |> str_replace(fixed("_"), "-") |> str_to_title()
      )
  }

  res |>
    relocate(ID, Index, Rarefaction, Diversity) |>
    update_provenance(se, list(analysis = "alpha diversity"))
}

test_alpha_diversity <- function(alpha_diversity, variable = "Group", two_sample_test = "wilcox", p_adjust_method = "fdr") {
  if (nrow(alpha_diversity) == 0L) {
    return(invisible())
  }

  alpha_diversity <- alpha_diversity |> mutate(across(any_of(variable), fortify))

  vi <-
    alpha_diversity |>
    dplyr::filter(Index == first(Index), Rarefaction == first(Rarefaction)) |>
    tibble_variable_info(variable)

  if (!vi[["testable"]]) {
    cli::cli_alert_warning(
      "{.field {provenance_as_short_title(alpha_diversity)}}: test skipped ({.var {variable}} must have at least two groups with at least two replicates each)"
    )
    return(invisible())
  }

  cli::cli_alert("{.field {provenance_as_short_title(alpha_diversity)}}: testing {.var {variable}}")

  results <-
    alpha_diversity |>
    dplyr::group_by(Index, Rarefaction) |>
    dplyr::group_map(\(x, y) bind_cols(y, pairwise_test(x, variable, "Diversity", two_sample_test, p_adjust_method)))

  results |>
    update_provenance(
      alpha_diversity,
      list(
        test = two_sample_test |> str_to_title() |> str_c(" Test"),
        `p-value correction` = p_adjust_method |> str_to_title() |> str_replace("^Fdr$", "FDR"),
        `variable of interest` = variable
      )
    )
}

format_alpha_diversity_test <- function(alpha_diversity_test_raw) {
  if (is.null(alpha_diversity_test_raw)) {
    return(invisible())
  }

  alpha_diversity_test_raw |>
    provenance_as_tibble() |>
    add_column(.y. = "Diversity") |>
    inner_join(dplyr::bind_rows(alpha_diversity_test_raw), by = ".y.") |>
    dplyr::select(!c(analysis, .y., p.adj.signif)) |>
    tibble::add_column(NAs = NA_character_, error = NA_character_) |>
    dplyr::rename(metric = Index, rarefaction = Rarefaction, `p-value` = p.adj, `uncorrected p-value` = p) |>
    dplyr::relocate(rarefaction, metric, .before = "variable of interest") |>
    dplyr::relocate(group1, group2, NAs, test, `p-value correction`, `p-value`, .after = "variable of interest") |>
    dplyr::relocate(`uncorrected p-value`, .after = last_col()) |>
    update_provenance(alpha_diversity_test_raw)
}

plot_alpha_diversity <- function(alpha_diversity, alpha_diversity_test_raw, variable_of_interest = "Group", theme = ggplot_theme()) {
  if (nrow(alpha_diversity) == 0L) {
    return(invisible())
  }

  alpha_diversity <- alpha_diversity |> mutate(across(any_of(variable_of_interest), fortify))

  vi <-
    alpha_diversity |>
    dplyr::filter(Index == first(Index), Rarefaction == first(Rarefaction)) |>
    tibble_variable_info(variable_of_interest)

  if (!vi[["testable"]]) {
    cli::cli_alert_warning(str_c(
      "{.field {provenance_as_short_title(alpha_diversity)}}: ",
      "plot skipped ({.var {variable_of_interest}} must have at least two groups with at least two replicates each)"
    ))
    return(invisible())
  }

  plot <- ggplot(
    data = alpha_diversity,
    mapping = aes(
      x = .data[[variable_of_interest]],
      y = Diversity
    )
  ) +
    geom_boxplot() +
    facet_grid(
      rows = vars(Index),
      cols = vars(Rarefaction),
      scales = "free_y"
    ) +
    labs(
      y = NULL
    ) +
    theme +
    theme(
      axis.text.x = element_text(angle = -90L, hjust = 0L, vjust = 0.5)
    )

  # assertion that plot and test belong together
  stopifnot(identical(
    alpha_diversity |> get_provenance() |> purrr::list_assign(`variable of interest` = variable_of_interest),
    alpha_diversity_test_raw |> get_provenance() |> purrr::list_assign(test = zap(), `p-value correction` = zap())
  ))

  # add p-values to the plot:
  if (!vec_is_empty(first(alpha_diversity_test_raw))) {
    ns_count <-
      alpha_diversity_test_raw |>
      map(\(x) {
        x |>
          dplyr::filter(p.adj.signif == "ns") |>
          vec_size()
      }) |>
      reduce(`+`)

    pvalue_data <-
      alpha_diversity_test_raw |>
      map(\(x) {
        attr(x, "args") <- list(
          data = x |> select(Index, Rarefaction) |> distinct() |> inner_join(alpha_diversity, by = join_by(Index, Rarefaction)),
          formula = as.formula(str_c("Diversity ~ ", variable_of_interest))
        )

        x |>
          rstatix::add_xy_position(x = variable_of_interest, step.increase = 0.4, scales = "free_y") |>
          dplyr::filter(p.adj.signif != "ns")
      }) |>
      dplyr::bind_rows()

    test <-
      alpha_diversity_test_raw |>
      get_provenance() |>
      chuck("test")
    correction <-
      alpha_diversity_test_raw |>
      get_provenance() |>
      chuck("p-value correction")

    plot <-
      plot +
      ggpubr::stat_pvalue_manual(
        pvalue_data,
        color = "blue",
        tip.length = 0L
      ) +
      labs(
        caption = pluralize(
          "{test}, {correction}-adjusted p-values: 0 **** 0.0001 *** 0.001 ** 0.01 * 0.05. {ns_count} non-significant comparison{?s} not shown."
        )
      ) +
      theme(
        plot.caption = element_text(color = "blue")
      )
  }

  plot |>
    update_provenance(
      alpha_diversity,
      list(
        aesthetics = list(by = variable_of_interest)
      )
    ) |>
    plot_titles(
      title = "alpha diversity analysis",
      analysis = zap()
    )
}
