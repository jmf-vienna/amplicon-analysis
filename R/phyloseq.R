as_phyloseq <- function(se) {
  ps <-
    se |>
    mia::makePhyloseqFromTreeSummarizedExperiment() |>
    microViz::tax_fix() |>
    microViz::phyloseq_validate()
  attr(ps, "provenance") <- S4Vectors::metadata(se)[["provenance"]]
  ps
}
