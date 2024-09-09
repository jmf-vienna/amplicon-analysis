make_summary_report <- function(provenance, pipeline_version, stats) {
  title <-
    provenance |>
    unlist() |>
    stringr::str_flatten(" ")
  glue::glue(
    "# {title} summary",
    "",
    "JMF downstream [amplicon analysis](https://github.com/jmf-vienna/amplicon-analysis) pipeline (beta) v{pipeline_version}",
    "",
    .sep = "\n",
    .trim = FALSE
  )
}
