make_col_data <- function(libraries) {
  libraries
}

make_row_data <- function(taxonomy) {
  taxonomy
}

make_assay_data <- function(counts) {
  names <- counts |> names()
  feature_column <- names[[1L]]
  library_column <- names[[2L]]
  count_column <- names[[3L]]

  counts |>
    tidyr::pivot_wider(
      id_cols = all_of(feature_column),
      names_from = all_of(library_column),
      names_sort = TRUE,
      values_from = all_of(count_column),
      values_fill = 0L
    ) |>
    tibble::column_to_rownames(feature_column) |>
    as.matrix()
}
