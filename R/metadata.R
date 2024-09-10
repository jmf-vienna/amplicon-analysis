variable_info <- function(x) {
  counts <- vctrs::vec_count(x, sort = "none")

  list(
    .length = length(x),
    .length_levels = nrow(counts),
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
}

ps_variable_info <- function(ps, variable_name) {
  data <- ps |> phyloseq::sample_data()

  if (!is.null(variable_name) && data |> rlang::has_name(variable_name)) {
    data |>
      dplyr::pull({{ variable_name }}) |>
      variable_info()
  } else {
    variable_info(list())
  }
}
