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
  removed_ids <- setdiff(rownames(before), rownames(after))
  removed <-
    SummarizedExperiment::rowData(before)[removed_ids, by, drop = FALSE] |>
    as_tibble("Feature_ID") |>
    dplyr::mutate(across(any_of("decontam_p_value"), \(x) cut(x, 0L:10L * 0.1, right = FALSE))) |>
    dplyr::mutate(across(where(is.numeric), as.character))

  removed_values <- removed |>
    dplyr::pull() |>
    jmf::uniques()
  cli::cli_alert("removed {nrow(removed)} feature{?s} via {by} filter{qty(removed_values)}{?/: /: }{removed_values}")

  new <- list(removed) |> set_names(by)
  S4Vectors::metadata(after) <- S4Vectors::metadata(after) |> list_modify(filtered_features = new)

  after
}

filtered_features_table <- function(se) {
  loadNamespace(class(se))

  se |>
    S4Vectors::metadata() |>
    purrr::pluck("filtered_features", .default = list()) |>
    map(\(x) {
      x |>
        dplyr::rename_with(\(x) "value", last_col()) |>
        dplyr::group_by(value) |>
        dplyr::summarise(count = dplyr::n(), features = Feature_ID |> str_flatten(" "))
    }) |>
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
