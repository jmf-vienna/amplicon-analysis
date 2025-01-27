make_col_data <- function(x) {
  suppressMessages(
    purrr::reduce(x, dplyr::inner_join)
  )
}

detect_taxonomy_ranks <- function(taxonomy) {
  ranks <- taxonomy |> names()
  cli::cli_alert("taxonomy ranks detected: {.val {ranks}}")
  ranks
}

make_row_data <- function(taxonomy, features_info) {
  dplyr::left_join(
    taxonomy,
    features_info,
    by = features_info |> first_name()
  )
}

make_assay_data <- function(counts) {
  counts <- counts
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

make_se <- function(counts, col_data, row_data, ranks, provenance) {
  sample_id_var <- col_data |> first_id_name()
  feature_id_var <- ranks |> dplyr::last()

  cli::cli_alert("ID variable names: sample = {.field {sample_id_var}} and feature = {.field {feature_id_var}}")

  col_data <-
    col_data |>
    dplyr::filter(.data[[sample_id_var]] %in% colnames(counts))

  row_data <-
    row_data |>
    dplyr::filter(.data[[feature_id_var]] %in% rownames(counts))

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
    rowData = row_data,
    metadata = list(
      taxonomy_ranks = ranks
    )
  ) |>
    set_provenance(provenance)
}

add_decontam <- function(se, negative_controls) {
  if (ncol(se) > 1L) {
    assay <- se |> SummarizedExperiment::assay()
    nc <- colnames(assay) %in% negative_controls

    decontam_result <- decontam::isContaminant(t(assay), neg = nc)
    stopifnot(identical(rownames(assay), rownames(decontam_result)))

    p <- decontam_result[["p"]]
  } else {
    cli::cli_alert_warning("skipped decontam (less than 2 libraries)")
    p <- NA_real_
  }

  SummarizedExperiment::rowData(se)[["decontam_p_value"]] <- p
  se
}

merge_cols <- function(se, by, keep_names, provenance = list()) {
  se <-
    mia::agglomerateByVariable(se, "cols", by) |>
    update_provenance(se, provenance)
  SummarizedExperiment::colData(se) <- SummarizedExperiment::colData(se)[keep_names]
  se
}

#### Overrides ####
# override these functions for project-specific logic

tidy_counts <- function(x) {
  x
}

tidy_samples <- function(x) {
  x
}

tidy_libraries <- function(x) {
  x
}

tidy_libraries_summary <- function(x) {
  x
}

tidy_taxonomy <- function(x) {
  x
}

tidy_features_info <- function(x) {
  x
}
