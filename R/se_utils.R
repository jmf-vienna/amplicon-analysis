tidy <- function(se) {
  se |> trim_empty()
}

trim_empty <- function(x, verbose = TRUE) {
  rs <- rowSums(SummarizedExperiment::assay(x))
  cs <- colSums(SummarizedExperiment::assay(x))

  if (any(rs == 0L)) {
    remove <- names(rs[rs == 0L])
    x <- x[!rownames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trimmed {length(remove)} empty feature(s)")
    }
  }

  if (any(cs == 0L)) {
    remove <- names(cs[cs == 0L])
    x <- x[, !colnames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trimmed {length(remove)} empty sample(s)")
    }
  }

  x
}

sample_id_var_name <- function(se) {
  se |>
    SummarizedExperiment::colData() |>
    names() |>
    head(1L)
}

feature_id_var_name <- function(se) {
  se |>
    SummarizedExperiment::rowData() |>
    names() |>
    tail(1L)
}

summary_as_row <- function(se) {
  summary <- se |> mia::summary()
  se |>
    provenance_as_tibble() |>
    tibble::add_column(
      sample = se |> sample_id_var_name(),
      feature = se |> feature_id_var_name(),
      samples = se |> SummarizedExperiment::colData() |> nrow(),
      features = summary[["features"]] |> dplyr::pull(total),
      total_counts = summary[["samples"]] |> dplyr::pull(total_counts),
      min_sample_counts = summary[["samples"]] |> dplyr::pull(min_counts),
      max_sample_counts = summary[["samples"]] |> dplyr::pull(max_counts),
      median_sample_counts = summary[["samples"]] |> dplyr::pull(median_counts)
    )
}

write_flattened <- function(se, file, assay_name = "counts") {
  assay_matrix <-
    se |>
    SummarizedExperiment::assay(assay_name) |>
    as.matrix()

  if (any(rowSums(assay_matrix) == 0L)) {
    empty <- rownames(assay_matrix)[rowSums(assay_matrix) == 0L]
    cli::cli_warn("rows with only zeros found in features {empty}")
  }
  if (any(colSums(assay_matrix) == 0L)) {
    empty <- colnames(assay_matrix)[colSums(assay_matrix) == 0L]
    cli::cli_warn("columns with only zeros found in samples {empty}")
  }

  assay_data <-
    assay_matrix |>
    as_tibble("ID") |>
    dplyr::mutate(across(where(is.numeric), as.character))

  row_data <-
    se |>
    SummarizedExperiment::rowData() |>
    as_tibble("ID")

  col_data <-
    se |>
    SummarizedExperiment::colData()
  col_data <-
    col_data |>
    purrr::map(as.character) |>
    as.data.frame(row.names = rownames(col_data)) |>
    as.matrix() |>
    t() |>
    as_tibble("ID")

  stopifnot(identical(
    assay_data |> dplyr::select(ID),
    row_data |> dplyr::select(ID)
  ))
  res <- dplyr::inner_join(assay_data, row_data, by = "ID")

  stopifnot(identical(
    assay_data |> names(),
    col_data |> names()
  ))
  res <- dplyr::bind_rows(col_data, res)

  cat(stringr::str_c("# ", se |> get_provenance() |> as_title(), "\n"), file = file)
  readr::write_tsv(res, file, na = "", append = TRUE, col_names = TRUE)

  invisible(file)
}

export_flattened <- function(se, dir_name, assay_name = "counts") {
  file <- fs::path(dir_name, se |> update_provenance(new = list(export = "flattened")) |> provenance_as_file_name(), ext = "tsv")
  prepare_export(file)

  cli::cli_alert("flattened data saved to {.file {file}}")
  write_flattened(se, file, assay_name)
}
