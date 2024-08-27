as_phyloseq <- function(se) {
  se |>
    trim_empty() |>
    mia::makePhyloseqFromTreeSummarizedExperiment() |>
    microViz::tax_fix() |>
    microViz::phyloseq_validate() |>
    update_provenance(se)
}
