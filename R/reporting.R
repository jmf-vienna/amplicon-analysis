make_summary_report <- function(provenance, pipeline_version, input_files, settings, config, stats) {
  title <-
    provenance |>
    unlist() |>
    stringr::str_flatten(" ")

  r_version <- stringr::str_c(R.Version()[["major"]], ".", R.Version()[["minor"]])

  packages <- c(
    "base",
    "targets",
    "tidyverse",
    "rstatix",
    "ggpubr",
    "SummarizedExperiment",
    "SingleCellExperiment",
    "mia",
    "vegan",
    "scuttle",
    "decontam",
    "phyloseq",
    "microbiome",
    "microViz"
  )

  citations <- stringr::str_c(
    "* `",
    packages,
    "` ",
    "v",
    packages |> purrr::map(packageVersion) |> purrr::map(as.character) |> unlist(),
    ": ",
    packages |> purrr::map(citation_text) |> unlist(),
    collapse = "\n"
  )

  input_files <- input_files |> discard(vec_is_empty)
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
    ifelse(
      isTRUE(all.equal(settings[["prevalence_filter"]], 0.0)) && isTRUE(all.equal(settings[["ra_filter"]], 0.0)),
      "",
      "Only features with a relative abundance of __more than {settings$ra_filter * 100} %__ in __more than {settings$prevalence_filter * 100} %__ of samples (excluding negative controls) were kept." # nolint
    ),
    ifelse(settings[["goods_coverage_min"]] > 0.0, "Minimum required Good’s Coverage was __{settings$goods_coverage_min}__.", ""),
    ifelse(settings[["yield_min"]] > 0L, "All samples with __less than {settings$yield_min} counts__ were removed.", ""),
    ifelse(settings[["yield_max"]] < Inf, "All samples with __more than {settings$yield_max} counts__ were removed.", ""),
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
    "Variable names:",
    "",
    "```yaml",
    "{yaml::as.yaml(settings$vars)}```",
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
