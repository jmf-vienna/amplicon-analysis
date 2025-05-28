add_alpha_diversity <- function(se, alpha_diversity_indexes, threshold, rarefaction_rounds = 10L) {
  min <- min_col_sum(se)
  if (min >= threshold) {
    se <-
      se |>
      mia::addAlpha(
        index = alpha_diversity_indexes,
        name = str_c(".alpha_diversity_", alpha_diversity_indexes)
      ) |>
      mia::addAlpha(
        index = alpha_diversity_indexes,
        name = str_c(".alpha_diversity_at_", min, "_", alpha_diversity_indexes),
        sample = min,
        niter = rarefaction_rounds
      )
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

plot_alpha_diversity <- function(alpha_diversity, variable, theme) {
  if (nrow(alpha_diversity) == 0L) {
    return(invisible())
  }

  alpha_diversity <- alpha_diversity |> mutate(across(any_of(variable), as.factor))

  vi <-
    alpha_diversity |>
    dplyr::filter(Index == first(Index), Rarefaction == first(Rarefaction)) |>
    tibble_variable_info(variable)

  if (!vi[["testable"]]) {
    cli::cli_alert_warning(
      "{.field {provenance_as_short_title(alpha_diversity)}}: plot skipped ({.var {variable}} must have at least two groups with at least two replicates each)"
    )
    return(invisible())
  }

  plot <- ggplot(
    data = alpha_diversity,
    mapping = aes(
      x = .data[[variable]],
      y = Diversity
    )
  ) +
    geom_boxplot() +
    facet_grid(
      rows = vars(Index),
      cols = vars(Rarefaction),
      scales = "free_y"
    ) +
    expand_limits(
      y = 0L
    ) +
    labs(
      y = NULL
    ) +
    theme +
    theme(
      axis.text.x = element_text(angle = -90L, hjust = 0L, vjust = 0.5)
    )

  plot |>
    update_provenance(alpha_diversity, list(
      aesthetics = list(by = variable)
    )) |>
    plot_titles(
      title = "alpha diversity analysis"
    )
}
