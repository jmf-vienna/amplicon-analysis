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
      cli::cli_alert("{.field {provenance_as_short_title(x)}}: trimmed {length(remove)} empty feature{?s}: {remove}")
    }
  }

  if (any(cs == 0L)) {
    remove <- names(cs[cs == 0L])
    x <- x[, !colnames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("{.field {provenance_as_short_title(x)}}: trimmed {length(remove)} empty sample{?s}: {remove}")
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
    as_full_tibble() |>
    dplyr::rename(sample_id = rowname, count = x)
}

min_col_sum <- function(se) {
  if (ncol(se) == 0L) {
    return(0L)
  }

  se |>
    col_sums() |>
    pull(count) |>
    min()
}

make_metrics <- function(se) {
  loadNamespace(class(se))

  features <-
    (SummarizedExperiment::assay(se) > 0L) |>
    colSums() |>
    as_full_tibble() |>
    dplyr::rename(sample_id = rowname, features = x)

  se |>
    provenance_as_tibble() |>
    tibble::add_column(
      tool = "JADA4",
      .after = "gene"
    ) |>
    dplyr::bind_cols(
      col_sums(se)
    ) |>
    inner_join(features, by = "sample_id") |>
    update_provenance(se, list(summary = "metrics"))
}

make_library_metrics <- function(se, library_id_var) {
  se |>
    make_metrics() |>
    dplyr::rename("{library_id_var}" := sample_id)
}

make_biosample_metrics <- function(se, biosample_id_var) {
  loadNamespace(class(se))

  n_libs <-
    se |>
    SummarizedExperiment::colData() |>
    as_full_tibble() |>
    dplyr::select(sample_id = 1L, libraries = .Number_of_libraries)

  se |>
    make_metrics() |>
    dplyr::inner_join(n_libs, by = "sample_id") |>
    dplyr::relocate(count, features, .after = last_col()) |>
    dplyr::rename("{biosample_id_var}" := sample_id)
}

make_metrics_summary <- function(library_metrics, biosample_metrics, library_id_var, biosample_id_var, se_summary) {
  dplyr::bind_rows(
    biosample_metrics |> dplyr::filter(resolution == "samples"),
    library_metrics
  ) |>
    dplyr::filter(count > 0L) |>
    dplyr::group_by(dplyr::across(!JMF_sample_ID:last_col())) |>
    dplyr::summarise(
      libraries = sum(libraries, na.rm = TRUE) + dplyr::n_distinct(.data[[library_id_var]], na.rm = TRUE),
      biosamples = dplyr::n_distinct(.data[[biosample_id_var]]),
      total_counts = sum(count),
      min_sample_counts = min(count),
      max_sample_counts = max(count),
      median_sample_counts = median(count),
      .groups = "drop"
    ) |>
    dplyr::left_join(se_summary) |>
    suppressMessages() |>
    dplyr::arrange(dplyr::across(!c(project:phase, libraries:last_col()), \(x) !is.na(x)))
}

summary_as_row <- function(se) {
  loadNamespace(class(se))

  sequence_length <- SummarizedExperiment::rowData(se)[["sequence_length"]]
  if (vec_is_empty(sequence_length)) {
    sequence_length <- NA_integer_
  }

  se |>
    provenance_as_tibble() |>
    tibble::add_column(
      features = nrow(se),
      min_sequence_length = min(sequence_length),
      max_sequence_length = max(sequence_length),
      median_sequence_length = median(sequence_length)
    )
}

plot_metrics <- function(data, hline_at, theme) {
  hline_at <- hline_at[hline_at > 0L & hline_at < Inf]
  id <- data |> first_id_name()

  plot <-
    ggplot(
      data = data,
      mapping = aes(
        x = .data[[id]] |> fct_rev(),
        y = count,
        label = str_c(scales::number(count), " ")
      )
    ) +
    geom_col() +
    geom_text(
      color = "white",
      hjust = 1L,
      size = 3L,
      family = font_family()
    ) +
    geom_hline(
      yintercept = hline_at,
      linetype = "dotted"
    ) +
    scale_y_log10(
      breaks = 10L^(0L:10L),
      minor_breaks = NULL,
      labels = scales::label_number()
    ) +
    coord_flip() +
    labs(
      x = id
    ) +
    theme

  plot |>
    update_provenance(
      data,
      list(
        aesthetics = "bar chart"
      )
    ) |>
    plot_titles(
      title = "yield",
      summary = zap()
    )
}

write_flattened <- function(se, file, assay_name = "counts") {
  loadNamespace(class(se))

  assay_matrix <-
    se |>
    SummarizedExperiment::assay(assay_name) |>
    as.matrix()

  if (any(rowSums(assay_matrix != 0L) == 0L)) {
    empty <- rownames(assay_matrix)[rowSums(assay_matrix != 0L) == 0L]
    cli::cli_warn("rows with only zeros found in features {empty}")
  }
  if (any(colSums(assay_matrix != 0L) == 0L)) {
    empty <- colnames(assay_matrix)[colSums(assay_matrix != 0L) == 0L]
    cli::cli_warn("columns with only zeros found in samples {empty}")
  }

  assay_data <-
    assay_matrix |>
    as_full_tibble("ID") |>
    dplyr::mutate(across(where(is.numeric), as.character))

  row_sums <-
    assay_matrix |>
    rowSums() |>
    as_full_tibble("ID") |>
    dplyr::rename(.row_sums = x)

  prevalence <-
    (assay_matrix > 0L) |>
    rowSums() |>
    as_full_tibble("ID") |>
    dplyr::rename(.prevalence = x)

  row_data <-
    se |>
    SummarizedExperiment::rowData() |>
    as_full_tibble("ID") |>
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

  features <-
    (assay_matrix != 0L) |>
    colSums() |>
    as.matrix() |>
    t() |>
    tibble::as_tibble() |>
    dplyr::mutate(across(everything(), as.character)) |>
    tibble::add_column(ID = ".Number_of_features", .before = 1L)

  col_data <-
    se |>
    SummarizedExperiment::colData()
  col_data <-
    col_data |>
    purrr::map(as.character) |>
    as.data.frame(row.names = rownames(col_data)) |>
    as.matrix() |>
    t() |>
    as_full_tibble("ID") |>
    dplyr::bind_rows(features) |>
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
  file <- fs::path(
    dir_name,
    stringr::str_c(
      provenance_as_file_name(se),
      "flattened",
      if (!identical(assay_name, "counts")) assay_name,
      sep = "_"
    ),
    ext = "tsv"
  )
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
