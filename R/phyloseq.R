as_phyloseq <- function(se) {
  se |>
    mia::makePhyloseqFromTreeSummarizedExperiment() |>
    microViz::tax_fix() |>
    microViz::phyloseq_validate()
}
