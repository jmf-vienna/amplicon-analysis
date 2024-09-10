keep_desirable_features <- function(se, config) {
  loadNamespace(class(se))

  purrr::reduce2(names(config), config, \(se, rank, values) {
    se_new <- se[SummarizedExperiment::rowData(se)[[rank]] %in% values]
    se_new <- filtered_features_helper(se, se_new, rank)
    se_new
  }, .init = se)
}

filter_undesirable_features <- function(se, config) {
  loadNamespace(class(se))

  se <- purrr::reduce2(names(config), config, \(se, rank, values) {
    se_new <- se[!SummarizedExperiment::rowData(se)[[rank]] %in% values]
    se_new <- filtered_features_helper(se, se_new, rank)
    se_new
  }, .init = se)
  se
}

filtered_features_helper <- function(before, after, rank) {
  features_before <- SummarizedExperiment::rowData(before)[[rank]]
  features_after <- SummarizedExperiment::rowData(after)[[rank]]

  n_removed <- length(features_before) - length(features_after)
  removed_features <- setdiff(features_before, features_after) |> jmf::uniques()
  cli::cli_alert("removed {n_removed} feature{?s} via {rank} filter")

  removed_features_count <-
    features_before[features_before %in% removed_features] |>
    vctrs::vec_count() |>
    tibble::as_tibble() |>
    dplyr::rename(feature = key, n = count)

  stopifnot(identical(n_removed, removed_features_count |> dplyr::pull(n) |> sum()))

  filtered_features <- list(
    list(
      n = n_removed,
      features = removed_features_count
    )
  )
  names(filtered_features) <- rank

  S4Vectors::metadata(after) <- S4Vectors::metadata(after) |> modifyList(list(filtered_features = filtered_features))

  after
}

filtered_features_table <- function(se) {
  loadNamespace(class(se))

  se |>
    S4Vectors::metadata() |>
    purrr::pluck("filtered_features", .default = list()) |>
    purrr::list_transpose() |>
    purrr::pluck("features", .default = tibble::tibble(feature = character(), n = integer())) |>
    dplyr::bind_rows(.id = "rank")
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
