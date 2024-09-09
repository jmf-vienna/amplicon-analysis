make_summary_report <- function(provenance, pipeline_version, stats) {
  title <-
    provenance |>
    unlist() |>
    stringr::str_flatten(" ")

  packages <- c(
    "targets", "tidyverse",
    "SummarizedExperiment", "SingleCellExperiment", "mia",
    "vegan", "scuttle",
    "phyloseq", "microViz"
  )

  citations <-
    packages |>
    purrr::map(citation_text) |>
    stringr::str_c("* `", packages, "`: ", text = _, collapse = "\n")

  glue::glue(
    "# {title} summary",
    "",
    "## Software versions",
    "",
    "* JMF downstream [amplicon analysis](https://github.com/jmf-vienna/amplicon-analysis) pipeline (beta) v{pipeline_version}",
    "* [`mia`](https://microbiome.github.io/mia/) v{packageVersion('mia')}",
    "* [`microViz`](https://david-barnett.github.io/microViz/) v{packageVersion('microViz')}",
    "",
    "## Citations",
    "",
    "{citations}",
    "",
    .sep = "\n",
    .trim = FALSE
  )
}

citation_text <- function(x) {
  withr::with_options(
    list(width = 10000),
    capture.output(print(citation(x), style = "text"))
  ) |>
    stringr::str_replace("\\*([0-9]+)\\*", "__\\1__")
}
