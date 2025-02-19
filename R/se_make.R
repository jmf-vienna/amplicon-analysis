make_col_data <- function(tables) {
  tables <- keep(tables, \(x) !vec_is_empty(x))
  suppressMessages(
    reduce(tables, dplyr::inner_join)
  )
}

taxonomy_fallback <- function(taxonomy, features_info) {
  if (ncol(taxonomy) == 0L) {
    cli::cli_alert_warning("falling back to empty taxonomy")
    fallback <-
      features_info |>
      dplyr::select(1L) |>
      tibble::add_column(
        Domain = "Unclassified",
        Phylum = NA_character_,
        Class = NA_character_,
        Order = NA_character_,
        Family = NA_character_,
        Genus = NA_character_,
        Species = NA_character_,
        .before = 1L
      )
    return(fallback)
  }

  taxonomy
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

get_failed_libraries <- function(se, negative_controls, pass_libraries_yield_min, failed_samples) {
  failed_libraries <-
    se |>
    col_sums() |>
    dplyr::filter(!ID %in% negative_controls, Sum < pass_libraries_yield_min) |>
    dplyr::pull(ID)

  if (!vec_is_empty(failed_libraries)) {
    cli::cli_alert("failed librar{?y/ies} (not enough yield): {.field {failed_libraries}}")
  }

  failed_libraries_via_sample <-
    se |>
    SummarizedExperiment::colData() |>
    as_tibble() |>
    dplyr::filter(.data[[biosample_id_var_name(se)]] %in% failed_samples) |>
    dplyr::pull(1L)

  if (!vec_is_empty(failed_libraries_via_sample)) {
    cli::cli_alert("failed librar{?y/ies} (marked via failed samples): {.field {failed_libraries_via_sample}}")
  }

  union(failed_libraries, failed_libraries_via_sample)
}

add_decontam <- function(se, negative_controls, failed_libraries = character()) {
  if (ncol(se) > 1L) {
    se_ready <- se[, !colnames(se) %in% failed_libraries]

    assay <- se_ready |> SummarizedExperiment::assay()
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
  loadNamespace("mia")

  lib_id <-
    se |>
    SummarizedExperiment::colData() |>
    first_id_name()

  lib_counts <-
    se |>
    SummarizedExperiment::colData() |>
    as_tibble() |>
    dplyr::group_by(across(all_of(by))) |>
    dplyr::summarise(
      .Number_of_libraries = dplyr::n(),
      .Libraries = .data[[lib_id]] |> str_flatten(" ")
    )

  se <-
    mia::agglomerateByVariable(se, "cols", by) |>
    update_provenance(se, provenance)

  SummarizedExperiment::colData(se) <-
    SummarizedExperiment::colData(se) |>
    as_tibble() |>
    dplyr::select(all_of(keep_names)) |>
    dplyr::inner_join(lib_counts, by = by) |>
    S4Vectors::DataFrame(row.names = colnames(se))

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

tidy_sublibraries <- function(x) {
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
