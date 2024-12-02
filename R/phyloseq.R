as_phyloseq <- function(se) {
  se |>
    mia::convertToPhyloseq() |>
    microViz::tax_fix() |>
    microViz::phyloseq_validate() |>
    update_provenance(se)
}
