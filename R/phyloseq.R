as_phyloseq <- function(se) {
  loadNamespace("mia")

  if (nrow(se) == 0L) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(se)}}: skipped because there are zero features")
    return(invisible())
  }

  se |>
    mia::convertToPhyloseq() |>
    microViz::tax_fix(anon_unique = FALSE, verbose = FALSE) |>
    microViz::phyloseq_validate() |>
    update_provenance(se)
}

export_ps <- function(ps, dir_name) {
  if (is.null(ps)) {
    return(invisible())
  }

  file <- fs::path(dir_name, ps |> update_provenance(new = list(export = "phyloseq")) |> provenance_as_file_name())
  c(
    write_rds(ps, fs::path(file, ext = "rds")),
    write_qs2(ps, fs::path(file, ext = "qs2"))
  )
}
