filter_contaminants <- function(se, decontam_threshold) {
  if (isTRUE(all.equal(decontam_threshold, 0.0))) {
    # keep everything
    return(se)
  }

  loadNamespace(class(se))

  p <- SummarizedExperiment::rowData(se)[["decontam_p_value"]]

  new_se <- se[is.na(p) | p >= decontam_threshold]
  filtered_features_helper(se, new_se, "decontam_p_value")
}

filter_by_length <- function(se, min = 0L, max = Inf) {
  if (identical(min, 0L) && identical(max, Inf)) {
    # keep everything
    return(se)
  }

  loadNamespace(class(se))

  l <- SummarizedExperiment::rowData(se)[["Sequence_length"]]

  new_se <- se[l >= min & l <= max]
  filtered_features_helper(se, new_se, "Sequence_length")
}

keep_desirable_features <- function(se, config) {
  if (vec_is_empty(config)) {
    # keep everything
    return(se)
  }

  loadNamespace(class(se))

  valid_ranks <- intersect(taxonomy_ranks(se), names(config))
  config <- config[valid_ranks]

  purrr::reduce2(names(config), config, \(se, rank, values) {
    se_new <- se[SummarizedExperiment::rowData(se)[[rank]] %in% values]
    se_new <- filtered_features_helper(se, se_new, rank)
    se_new
  }, .init = se)
}

filter_undesirable_features <- function(se, config) {
  if (vec_is_empty(config)) {
    # keep everything
    return(se)
  }

  loadNamespace(class(se))

  valid_ranks <- intersect(taxonomy_ranks(se), names(config))
  config <- config[valid_ranks]

  se <- purrr::reduce2(names(config), config, \(se, rank, values) {
    se_new <- se[!SummarizedExperiment::rowData(se)[[rank]] %in% values]
    se_new <- filtered_features_helper(se, se_new, rank)
    se_new
  }, .init = se)
  se
}

filtered_features_helper <- function(before, after, by) {
  features_before <- SummarizedExperiment::rowData(before)[[by]]
  features_after <- SummarizedExperiment::rowData(after)[[by]]

  if (by == "decontam_p_value") {
    features_before <- features_before |> cut(0L:10L * 0.1, right = FALSE)
    features_after <- features_after |> cut(0L:10L * 0.1, right = FALSE)
  } else if (by == "Sequence_length") {
    features_before <- features_before |> as.character()
    features_after <- features_after |> as.character()
  }

  n_removed <- length(features_before) - length(features_after)
  removed_features <- setdiff(features_before, features_after) |> jmf::uniques()
  cli::cli_alert("removed {n_removed} feature{?s} via {by} filter{qty(removed_features)}{?/: /: }{removed_features}")

  removed_features_count <-
    features_before[features_before %in% removed_features] |>
    vctrs::vec_count() |>
    tibble::as_tibble() |>
    dplyr::rename(value = key) |>
    dplyr::arrange(value)

  stopifnot(identical(n_removed, removed_features_count |> dplyr::pull(count) |> sum()))

  filtered_features <- list(
    list(
      n = n_removed,
      features = removed_features_count
    )
  )
  names(filtered_features) <- by

  S4Vectors::metadata(after) <- S4Vectors::metadata(after) |> modifyList(list(filtered_features = filtered_features))

  after
}

filtered_features_table <- function(se) {
  loadNamespace(class(se))

  se |>
    S4Vectors::metadata() |>
    purrr::pluck("filtered_features", .default = list()) |>
    purrr::list_transpose() |>
    purrr::pluck("features", .default = tibble::tibble(value = character(), count = integer())) |>
    dplyr::bind_rows(.id = "filter")
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
