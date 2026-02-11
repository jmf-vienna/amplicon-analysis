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

  l <- SummarizedExperiment::rowData(se)[["sequence_length"]]

  new_se <- se[l >= min & l <= max]
  filtered_features_helper(se, new_se, "sequence_length")
}

keep_desirable_features <- function(se, config) {
  if (vec_is_empty(config)) {
    # keep everything
    return(se)
  }

  loadNamespace(class(se))

  valid_ranks <- intersect(taxonomy_ranks(se), names(config))
  config <- config[valid_ranks]

  purrr::reduce2(
    names(config),
    config,
    \(se, rank, values) {
      se_new <- se[SummarizedExperiment::rowData(se)[[rank]] %in% values]
      se_new <- filtered_features_helper(se, se_new, rank)
      se_new
    },
    .init = se
  )
}

filter_undesirable_features <- function(se, config) {
  if (vec_is_empty(config)) {
    # keep everything
    return(se)
  }

  loadNamespace(class(se))

  valid_ranks <- intersect(taxonomy_ranks(se), names(config))
  config <- config[valid_ranks]

  se <- purrr::reduce2(
    names(config),
    config,
    \(se, rank, values) {
      se_new <- se[!SummarizedExperiment::rowData(se)[[rank]] %in% values]
      se_new <- filtered_features_helper(se, se_new, rank)
      se_new
    },
    .init = se
  )
  se
}

filtered_features_helper <- function(before, after, by) {
  removed_ids <- setdiff(rownames(before), rownames(after))
  removed <-
    SummarizedExperiment::rowData(before)[removed_ids, by, drop = FALSE] |>
    as_full_tibble("Feature_ID") |>
    dplyr::mutate(
      across(any_of("decontam_p_value"), \(x) cut(x, 0L:10L * 0.1, right = FALSE)),
      across(where(is.numeric), as.character)
    )

  removed_values <- removed |>
    dplyr::pull() |>
    jmf::uniques()
  cli::cli_alert("removed {nrow(removed)} feature{?s} via {by} filter{qty(removed_values)}{?/: /: }{removed_values}")

  new <- list(removed) |> set_names(by)
  S4Vectors::metadata(after) <- S4Vectors::metadata(after) |> list_modify(filtered_features = new)

  after
}

make_filtered_features_table <- function(se) {
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
    dplyr::bind_rows(.id = "filter") |>
    dplyr::bind_rows(
      # default empty table
      tibble::tibble(filter = NA_character_, value = NA_character_, count = NA_integer_, features = NA_character_, .rows = 0L)
    )
}

remove_cols <- function(se, cols) {
  loadNamespace(class(se))
  se[, !colnames(se) %in% cols]
}

filter_samples_by_sum <- function(se, min = 0L, max = Inf) {
  keep <-
    se |>
    col_sums() |>
    dplyr::filter(count >= min, count <= max) |>
    dplyr::pull(sample_id)

  res <- se[, keep]

  if (min > 0L) {
    res <- res |> update_provenance(new = list("sample filter" = list("≥" = min)))
  }

  if (max < Inf) {
    res <- res |> update_provenance(new = list("sample filter" = list("≤" = min)))
  }

  res
}

filter_prevalence_ra <- function(se, prev = 0.0, ra = 0.0, negative_controls = character()) {
  if (isTRUE(all.equal(prev, 0.0)) && isTRUE(all.equal(ra, 0.0))) {
    return(se)
  }

  loadNamespace(class(se))

  ready <- se[, !colnames(se) %in% negative_controls]

  n_prev <- floor(ncol(ready) * prev)

  above_ra <- rowSums(SummarizedExperiment::assay(ready, "relabundance") > ra)
  sel <- above_ra > n_prev

  new_se <-
    se[sel] |>
    update_provenance(new = list(state = "filtered"))

  filter_name <- str_c("n_samples(relative abundance > ", ra * 100L, " %)")
  SummarizedExperiment::rowData(se)[filter_name] <- above_ra

  filtered_features_helper(se, new_se, filter_name)
}

filter_samples_by_goods_coverage <- function(se, min = 0.0) {
  goods_cov <-
    se |>
    SummarizedExperiment::assay() |>
    get_goods_coverage()
  keep <- names(goods_cov[goods_cov >= min])

  res <- se[, keep]

  if (min > 0.0) {
    res <- res |> update_provenance(new = list("sample filter" = list("Good’s Coverage ≥" = min)))
  }

  res
}

make_filtered_samples_table <- function(se_pairs) {
  loadNamespace(class(se_pairs |> chuck(1L, 1L)))

  map(se_pairs, \(se_pair) {
    se1 <- se_pair |> chuck(1L)
    se2 <- se_pair |> chuck(2L)

    se2 |>
      provenance_as_tibble() |>
      tibble::add_column(
        removed = setdiff(colnames(se1), colnames(se2)) |> str_flatten(" ")
      )
  }) |>
    rev() |>
    dplyr::bind_rows() |>
    dplyr::arrange(desc(dplyr::row_number()))
}
