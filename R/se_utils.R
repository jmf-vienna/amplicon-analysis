trim_empty <- function(x, verbose = TRUE) {
  rs <- rowSums(SummarizedExperiment::assay(x))
  cs <- colSums(SummarizedExperiment::assay(x))

  if (any(rs == 0L)) {
    remove <- names(rs[rs == 0L])
    x <- x[!rownames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trim_empty: Trimmed {length(remove)} empty row(s).")
    }
  }

  if (any(cs == 0L)) {
    remove <- names(cs[cs == 0L])
    x <- x[, !colnames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trim_empty: Trimmed {length(remove)} empty column(s).")
    }
  }

  x
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
    as.data.frame() |>
    tibble::rownames_to_column("FeatureID") |>
    tibble::as_tibble() |>
    dplyr::mutate(across(where(is.numeric), as.character))

  row_data <-
    se |>
    SummarizedExperiment::rowData() |>
    as.data.frame() |>
    tibble::rownames_to_column("FeatureID") |>
    tibble::as_tibble()

  col_data <-
    se |>
    SummarizedExperiment::colData()
  col_data <-
    col_data |>
    purrr::map(as.character) |>
    as.data.frame(row.names = rownames(col_data)) |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("FeatureID") |>
    tibble::as_tibble()

  stopifnot(identical(
    assay_data |> dplyr::select(FeatureID),
    row_data |> dplyr::select(FeatureID)
  ))
  res <- dplyr::inner_join(assay_data, row_data, by = "FeatureID")

  stopifnot(identical(
    assay_data |> names(),
    col_data |> names()
  ))
  res <- dplyr::bind_rows(col_data, res)

  res <- res |> dplyr::rename(ID = FeatureID)

  cat(stringr::str_c("# ", se |> get_provenance() |> as_title(), "\n"), file = file)
  readr::write_tsv(res, file, na = "", append = TRUE, col_names = TRUE)

  invisible(file)
}

export_flattened <- function(se, dir_name, assay_name = "counts") {
  file <- fs::path(dir_name, se |> as_file_name(), ext = "tsv")
  prepare_export(file)

  write_flattened(se, file, assay_name)
}
