trim_empty <- function(x, verbose = TRUE) {
  rs <- rowSums(SummarizedExperiment::assay(x))
  cs <- colSums(SummarizedExperiment::assay(x))

  if (any(rs == 0L)) {
    remove <- names(rs[rs == 0L])
    x <- x[!rownames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trim_empty: Trimmed {length(remove)} empty row(s).")
    }
  }

  if (any(cs == 0L)) {
    remove <- names(cs[cs == 0L])
    x <- x[, !colnames(x) %in% remove]

    if (verbose) {
      cli::cli_alert("trim_empty: Trimmed {length(remove)} empty column(s).")
    }
  }

  x
}
