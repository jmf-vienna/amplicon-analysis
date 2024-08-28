keep_desirable_features <- function(se, config) {
  se |>
    SummarizedExperiment::subset(
      Domain %in% config[["Domain"]]
    )
}

filter_undesirable_features <- function(se, config) {
  se |>
    SummarizedExperiment::subset(!(
      Order %in% config[["Order"]] |
        Family %in% config[["Family"]] |
        ASV %in% config[["ASV"]]
    ))
}

filter_samples_by_sum <- function(se, min = 0L, max = Inf) {
  keep <-
    se |>
    scuttle::perCellQCMetrics(use_altexps = FALSE) |>
    as_tibble() |>
    dplyr::filter(sum >= min, sum <= max) |>
    dplyr::pull(1L)

  res <- se[, keep]

  if (min > 0L) {
    res <- res |> update_provenance(new = list("sample filter" = list("≥" = min)))
  }

  if (max < Inf) {
    res <- res |> update_provenance(new = list("sample filter" = list("≤" = min)))
  }

  res
}
