variable_info <- function(x, label) {
  exists <- !is.null(x)
  if (!exists) {
    x <- list()
  }

  counts <- vctrs::vec_count(x, sort = "none")
  n_na <- sum(is.na(x))

  # nolint start: nrow_subset_linter.
  list(
    .label = label,
    .length = length(x),
    .length_levels = nrow(counts),
    .NAs = n_na,
    exists = exists,
    has_na = n_na >= 1L,
    all_na = all(is.na(x)),
    # more than one value:
    multiple = length(x) > 1L,
    # more than one level:
    multiple_levels = nrow(counts) > 1L,
    # at least one group with at least two replicates each
    duplicates = counts |>
      dplyr::filter(count >= 2L) |>
      nrow() >=
      1L,
    # at least two groups with at least two replicates each
    testable = counts |>
      dplyr::filter(count >= 2L) |>
      nrow() >=
      2L,
    # all groups must have at least two replicates each (with at least two groups)
    pairwise_testable = counts |>
      dplyr::filter(count >= 2L) |>
      nrow() ==
      max(2L, nrow(counts))
  )
  # nolint end
}

tibble_variable_info <- function(data, variable_name) {
  if (!is.null(variable_name) && rlang::has_name(data, variable_name)) {
    data |>
      dplyr::pull({{ variable_name }}) |>
      variable_info(variable_name)
  } else {
    variable_info(NULL, variable_name)
  }
}

ps_variable_info <- function(ps, variable_name) {
  ps |>
    phyloseq::sample_data() |>
    tibble_variable_info(variable_name)
}
