as_phyloseq <- function(se) {
  se |>
    mia::convertToPhyloseq() |>
    microViz::tax_fix() |>
    microViz::phyloseq_validate() |>
    update_provenance(se)
}

export_ps <- function(se, dir_name) {
  file <- fs::path(dir_name, se |> update_provenance(new = list(export = "phyloseq")) |> provenance_as_file_name(), ext = "rds")
  write_rds(se, file)
}
