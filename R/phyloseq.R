as_phyloseq <- function(se, assay_type = "counts") {
  if (is_too_large(se)) {
    return()
  }

  if (!assay_type %in% SummarizedExperiment::assayNames(se)) {
    return()
  }

  loadNamespace("mia")

  if (nrow(se) == 0L) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(se)}}: skipped because there are zero features")
    return(invisible())
  }

  se |>
    mia::convertToPhyloseq(assay.type = assay_type) |>
    microViz::tax_fix(anon_unique = FALSE, verbose = FALSE) |>
    microViz::phyloseq_validate() |>
    update_provenance(se, list(assay_type = assay_type))
}

export_ps <- function(ps, dir_name) {
  if (is.null(ps)) {
    return(invisible())
  }

  file <- fs::path(dir_name, ps |> update_provenance(new = list(export = "phyloseq")) |> provenance_as_file_name())

  write_rd(ps, file)
}
