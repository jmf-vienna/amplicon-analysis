variable_info <- function(x) {
  exists <- !is.null(x)
  if (!exists) {
    x <- list()
  }

  counts <- vctrs::vec_count(x, sort = "none")
  n_na <- sum(is.na(x))

  # nolint start: nrow_subset_linter.
  list(
    .length = length(x),
    .length_levels = nrow(counts),
    .NAs = n_na,
    exists = exists,
    has_na = n_na > 1L,
    all_na = all(is.na(x)),
    # more than one value:
    multiple = length(x) > 1L,
    # more than one level:
    multiple_levels = nrow(counts) > 1L,
    # at least one group with at least two replicates each
    duplicates =
      counts |>
        dplyr::filter(count >= 2L) |>
        nrow()
      >= 1L,
    # at least two groups with at least two replicates each
    testable =
      counts |>
        dplyr::filter(count >= 2L) |>
        nrow()
      >= 2L
  )
  # nolint end
}

ps_variable_info <- function(ps, variable_name) {
  data <- ps |> phyloseq::sample_data()

  if (!is.null(variable_name) && data |> rlang::has_name(variable_name)) {
    data |>
      dplyr::pull({{ variable_name }}) |>
      variable_info()
  } else {
    variable_info(NULL)
  }
}
