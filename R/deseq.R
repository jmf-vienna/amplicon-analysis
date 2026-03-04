deseq <- function(se) {
  se |>
    purrr::pluck(attr_getter("analyze"), "category") |>
    purrr::map(\(x) run_deseq(se, x))
}

run_deseq <- function(se, var, pseudocount = 1L, min_features = 3L, alpha = 0.05) {
  loadNamespace(class(se))

  if (nrow(se) < min_features) {
    cli::cli_alert_warning("{.field {provenance_as_short_title(se)}}: DESeq2 skipped on {.var {var}}: ≥{min_features} feature{?s} required (found {nrow(se)})")
    return(invisible())
  }

  SummarizedExperiment::colData(se)[[var]] <-
    SummarizedExperiment::colData(se)[[var]] |>
    fortify()

  vi <- se_variable_info(se, var)

  if (!vi[["two_levels"]]) {
    cli::cli_alert_danger("{.field {provenance_as_short_title(se)}}: DESeq2 skipped on {.var {var}}: must have exactly two levels")
    return(invisible())
  }

  if (!vi[["testable"]]) {
    cli::cli_alert_warning(
      "{.field {provenance_as_short_title(se)}}: DESeq2 skipped on ({.var {var}}: must have at least two groups with at least two replicates each)"
    )
    return(invisible())
  }

  cli::cli_alert("{.field {provenance_as_short_title(se)}}: running DESeq2 with {.var {var}}")

  deseq <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(se) + pseudocount,
    colData = SummarizedExperiment::colData(se),
    design = as.formula(stringr::str_c("~ ", var))
  )

  deseq <- DESeq2::DESeq(deseq, quiet = !interactive())

  groups <- levels(SummarizedExperiment::colData(se)[[var]])

  row_data <-
    se |>
    SummarizedExperiment::rowData() |>
    as_full_tibble("Feature_ID") |>
    dplyr::select(Feature_ID, Lineage)

  res <-
    DESeq2::results(deseq, alpha = alpha) |>
    as_full_tibble("Feature_ID") |>
    dplyr::select(Feature_ID, `log2 fold change` = log2FoldChange, `p-value` = padj, `uncorrected p-value` = pvalue) |>
    tibble::add_column(`variable of interest` = var, group1 = groups[[1L]], group2 = groups[[2L]], .before = "Feature_ID") |>
    bind_cols(x = provenance_as_tibble(se), y = _) |>
    dplyr::left_join(row_data, by = "Feature_ID")
}

collect_deseq_results <- function(deseq_raw_results) {
  res <-
    deseq_raw_results |>
    bind_rows()

  if (vec_is_empty(res)) {
    res <- res |>
      dplyr::mutate(
        `log2 fold change` = NA_real_,
        `p-value` = NA_real_,
        `uncorrected p-value` = NA_real_
      )
  }

  res
}

filter_deseq_results <- function(deseq_combined_results, log2_fold_change = 2.0, p_value = 0.05) {
  deseq_combined_results |>
    dplyr::filter(abs(`log2 fold change`) >= log2_fold_change, `p-value` <= p_value) |>
    # `signif` for reproducibility. Keeping six digits is already overkill.
    dplyr::mutate(`log2 fold change` = signif(`log2 fold change`))
}

split_tibble_by_rank <- function(x) {
  x |>
    group_by(rank) |>
    dplyr::group_map(\(d, k) update_provenance(d, x, new = as.list(k)))
}

plot_deseq <- function(plot_data, theme) {
  if (vec_is_empty(plot_data)) {
    return(invisible())
  }

  n_features <-
    plot_data |>
    pull(Feature_ID) |>
    unique() |>
    vec_size()

  n_tests <-
    plot_data |>
    distinct(across(!Feature_ID:last_col())) |>
    nrow()

  lin_len <-
    plot_data |>
    pull(Lineage) |>
    unique() |>
    nchar() |>
    max()

  p <-
    ggplot(
      data = plot_data,
      mapping = aes(
        x = str_c(`variable of interest`, ":\n", group2, "\nvs\n", group1),
        y = Lineage |> fct_rev(),
        size = abs(`log2 fold change`),
        fill = `log2 fold change` |> sign() |> factor(c("-1", "1")) |> fct_recode("↓" = "-1", "↑" = "1")
      )
    ) +
    geom_point(
      shape = 21L
    ) +
    scale_fill_manual(
      values = c("↓" = "white", "↑" = "black")
    ) +
    facet_grid(
      cols = vars(subset |> str_replace(fixed(": "), ":\n")),
      scales = "free",
      space = "free"
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Change",
      size = "log\u2082FC"
    ) +
    theme +
    theme(
      legend.position = "bottom"
    )

  p <- p |>
    update_provenance(plot_data) |>
    plot_titles(title = "DESeq2 differential abundance test", test = zap())

  attr(p, "output") <- list(
    width = 0.06 * lin_len + 1.0 * max(n_tests, 1L),
    height = 2.5 + 0.2 * max(n_features, 1L)
  )

  p
}
