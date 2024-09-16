make_summary_report <- function(provenance, pipeline_version, input_files, settings, config, stats) {
  title <-
    provenance |>
    unlist() |>
    stringr::str_flatten(" ")

  r_version <- stringr::str_c(R.Version()[["major"]], ".", R.Version()[["minor"]])

  packages <- c(
    "base",
    "targets", "tidyverse",
    "SummarizedExperiment", "SingleCellExperiment", "mia",
    "vegan", "scuttle",
    "phyloseq", "microViz"
  )

  citations <-
    packages |>
    purrr::map(citation_text) |>
    stringr::str_c("* `", packages, "`: ", text = _, collapse = "\n")

  input_files_md <- stringr::str_c("* ", names(input_files), ": `", input_files, "`", collapse = "\n")

  glue::glue(
    "# {title} summary",
    "",
    "## Filter settings",
    "",
    "__Desirable__ taxa, i.e., keep only these:",
    "",
    "```yaml",
    "{yaml::as.yaml(settings$desirables)}```",
    "",
    "__Undesirable__ taxa, i.e., remove these:",
    "",
    "```yaml",
    "{yaml::as.yaml(settings$undesirables)}```",
    "",
    "All samples with __less than {settings$yield_min} counts__ were removed.",
    "",
    "## Software versions",
    "",
    "* JMF downstream [amplicon analysis](https://github.com/jmf-vienna/amplicon-analysis) pipeline v{pipeline_version}",
    "* [`mia`](https://microbiome.github.io/mia/) v{packageVersion('mia')}",
    "* [`microViz`](https://david-barnett.github.io/microViz/) v{packageVersion('microViz')}",
    "* [Bioconductor](https://www.bioconductor.org/) v{packageVersion('BiocVersion')}",
    "* R v{r_version}",
    "",
    "## Citations",
    "",
    "{citations}",
    "",
    "## Input files",
    "",
    "{input_files_md}",
    "",
    "## Complete config",
    "",
    "```yaml",
    "{yaml::as.yaml(config)}```",
    "",
    .sep = "\n",
    .trim = FALSE
  )
}

citation_text <- function(x) {
  withr::with_options(
    list(width = 10000L),
    capture.output(print(citation(x), style = "text"))
  ) |>
    stringr::str_replace("\\*([0-9]+)\\*", "__\\1__") # nolint: nonportable_path_linter.
}

make_previous_summary_rows <- function(x, provenance) {
  sample_id_var_name <-
    x |>
    names() |>
    head(1L)

  x |>
    dplyr::group_by(phase) |>
    dplyr::summarise(
      samples = dplyr::n_distinct(.data[[sample_id_var_name]]),
      total_counts = sum(count),
      min_sample_counts = min(count),
      max_sample_counts = max(count),
      median_sample_counts = median(count)
    ) |>
    tibble::add_column(
      project = provenance |> purrr::chuck("project"),
      gene = provenance |> purrr::chuck("gene"),
      resolution = "libraries",
      state = "crude",
      .before = "phase"
    ) |>
    tibble::add_column(
      sample = sample_id_var_name,
      .before = "samples"
    )
}
