make_col_data <- function(x) {
  suppressMessages(
    purrr::reduce(x, dplyr::inner_join)
  )
}

# override for custom user logic
tidy_taxonomy <- function(x) {
  x
}

make_row_data <- function(taxonomy) {
  taxonomy |>
    tidy_taxonomy()
}

make_assay_data <- function(counts) {
  names <- counts |> names()
  feature_column <- names[[1L]]
  samples_column <- names[[2L]]
  count_column <- names[[3L]]

  counts |>
    tidyr::pivot_wider(
      id_cols = all_of(feature_column),
      names_from = all_of(samples_column),
      names_sort = TRUE,
      values_from = all_of(count_column),
      values_fill = 0L
    ) |>
    tibble::column_to_rownames(feature_column) |>
    as.matrix()
}

make_se <- function(counts, col_data, row_data, provenance) {
  sample_id_var <-
    col_data |>
    names() |>
    head(1L)

  feature_id_var <-
    row_data |>
    names() |>
    tail(1L)

  col_data <-
    col_data |>
    dplyr::filter(.data[[sample_id_var]] %in% colnames(counts))

  stopifnot(
    identical(
      col_data |> dplyr::pull(sample_id_var),
      counts |> colnames()
    ),
    identical(
      row_data |> dplyr::pull(feature_id_var),
      counts |> rownames()
    )
  )

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = col_data,
    rowData = row_data
  ) |>
    set_provenance(provenance)
}

merge_cols <- function(se, by, keep_names, provenance = list()) {
  se <-
    mia::mergeCols(se, se[[by]]) |>
    update_provenance(se, provenance)
  SummarizedExperiment::colData(se) <- SummarizedExperiment::colData(se)[keep_names]
  se
}
