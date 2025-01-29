tidy <- function(se) {
  se |> trim_empty()
}

taxonomy_ranks <- function(se) {
  se |>
    S4Vectors::metadata() |>
    chuck("taxonomy_ranks")
}

trim_empty <- function(x, verbose = TRUE) {
  loadNamespace(class(x))

  rs <- rowSums(SummarizedExperiment::assay(x))
  cs <- colSums(SummarizedExperiment::assay(x))

  if (any(rs == 0L)) {
    remove <- names(rs[rs == 0L])
    x <- x[!rownames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trimmed {length(remove)} empty feature{?s}: {remove}")
    }
  }

  if (any(cs == 0L)) {
    remove <- names(cs[cs == 0L])
    x <- x[, !colnames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trimmed {length(remove)} empty sample{?s}: {remove}")
    }
  }

  x
}

sample_id_var_name <- function(se) {
  se |>
    SummarizedExperiment::colData() |>
    first_name()
}

library_id_var_name <- function(se) {
  se |>
    SummarizedExperiment::colData() |>
    names() |>
    stringr::str_subset("(lib|Lib|LIB|library|Library|LIBRARY).+ID$") |>
    dplyr::first()
}

biosample_id_var_name <- function(se) {
  se |>
    SummarizedExperiment::colData() |>
    names() |>
    stringr::str_subset("(sample|Sample|SAMPLE).+ID$") |>
    dplyr::first()
}

feature_id_var_name <- function(se) {
  se |>
    taxonomy_ranks() |>
    dplyr::last()
}

col_sums <- function(se) {
  loadNamespace(class(se))

  se |>
    SummarizedExperiment::assay() |>
    colSums() |>
    as_tibble() |>
    dplyr::rename(ID = rowname, Sum = x)
}

summary_as_row <- function(se) {
  loadNamespace(class(se))

  col_data <-
    se |>
    SummarizedExperiment::colData() |>
    as_tibble()

  biosample_id_var_name <- se |> biosample_id_var_name()
  n_biosamples <-
    col_data |>
    dplyr::pull(biosample_id_var_name) |>
    vec_unique_count()

  library_id_var_name <- se |> library_id_var_name()
  n_libraries <- if (col_data |> tibble::has_name(".Number_of_libraries")) {
    col_data |>
      dplyr::pull(.Number_of_libraries) |>
      sum()
  } else {
    col_data |>
      dplyr::pull(library_id_var_name) |>
      vec_unique_count()
  }

  sample_counts <-
    se |>
    col_sums() |>
    dplyr::pull(Sum)
  if (vec_is_empty(sample_counts)) sample_counts <- NA_integer_

  sequence_length <- SummarizedExperiment::rowData(se)[["Sequence_length"]]
  if (vec_is_empty(sequence_length)) sequence_length <- NA_integer_

  se |>
    provenance_as_tibble() |>
    tibble::add_column(
      tool = "JADA4",
      sample = se |> sample_id_var_name(),
      feature = se |> feature_id_var_name(),
      biosamples = n_biosamples,
      libraries = n_libraries,
      features = nrow(se),
      total_counts = sum(sample_counts),
      min_sample_counts = min(sample_counts),
      max_sample_counts = max(sample_counts),
      median_sample_counts = median(sample_counts),
      min_sequence_length = min(sequence_length),
      max_sequence_length = max(sequence_length),
      median_sequence_length = median(sequence_length)
    )
}

write_flattened <- function(se, file, assay_name = "counts") {
  loadNamespace(class(se))

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

  row_sums <-
    assay_matrix |>
    rowSums() |>
    as_tibble("ID") |>
    dplyr::rename(.row_sums = x)

  prevalence <-
    (assay_matrix > 0L) |>
    rowSums() |>
    as_tibble("ID") |>
    dplyr::rename(.prevalence = x)

  row_data <-
    se |>
    SummarizedExperiment::rowData() |>
    as_tibble("ID") |>
    dplyr::left_join(row_sums, by = "ID") |>
    dplyr::left_join(prevalence, by = "ID") |>
    dplyr::relocate(.row_sums:last_col(), .after = ID)

  col_sums <-
    assay_matrix |>
    colSums() |>
    as.matrix() |>
    t() |>
    tibble::as_tibble() |>
    dplyr::mutate(across(everything(), as.character)) |>
    tibble::add_column(ID = ".col_sums", .before = 1L)

  col_data <-
    se |>
    SummarizedExperiment::colData()
  col_data <-
    col_data |>
    purrr::map(as.character) |>
    as.data.frame(row.names = rownames(col_data)) |>
    as.matrix() |>
    t() |>
    as_tibble("ID") |>
    dplyr::bind_rows(col_sums)

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

export_se <- function(se, dir_name) {
  file <- fs::path(dir_name, se |> update_provenance(new = list(export = "SE")) |> provenance_as_file_name())
  c(
    write_rds(se, fs::path(file, ext = "rds")),
    write_qs2(se, fs::path(file, ext = "qs2"))
  )
}
