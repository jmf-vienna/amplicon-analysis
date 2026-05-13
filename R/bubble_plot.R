fill_unclassified <- function(se, value = "unclassified", species_value = "sp.", parentheses = TRUE) {
  loadNamespace(class(se))

  if (isTRUE(parentheses)) {
    prefix <- "("
    suffix <- ")"
  } else {
    prefix <- ""
    suffix <- ""
  }

  row_data <-
    se |>
    SummarizedExperiment::rowData() |>
    as_full_tibble()

  prev_rank <- NA
  for (rank in taxonomy_ranks(se)) {
    if (is.na(prev_rank)) {
      row_data <-
        row_data |>
        dplyr::mutate(
          "{rank}" := .data[[rank]] |> tidyr::replace_na(stringr::str_c(prefix, value, suffix))
        )
    } else {
      row_data <-
        row_data |>
        mutate(
          "{rank}" := dplyr::coalesce(
            .data[[rank]],
            if_else(
              stringr::str_detect(.data[[prev_rank]], value),
              .data[[prev_rank]],
              stringr::str_c(prefix, value, " ", .data[[prev_rank]], suffix)
            )
          )
        )
    }

    if (rlang::is_string(species_value) && prev_rank == "Genus" && rank == "Species") {
      row_data <-
        row_data |>
        mutate(
          Species = if_else(
            # if: Genus IS NOT unclassified & Species IS unclassified
            !stringr::str_detect(Genus, value) & stringr::str_detect(Species, value),
            # then: Genus sp.
            stringr::str_c(prefix, Genus, " ", species_value, suffix),
            # else: leave as is
            Species
          )
        )
    }

    prev_rank <- rank
  }

  stopifnot(identical(row_data[["rowname"]], rownames(se)))
  SummarizedExperiment::rowData(se) <- tibble::column_to_rownames(row_data)

  se
}

smart_agglomerate <- function(se, min_abundance = 0.01, min_prevalence = 2L, remove_zeros = TRUE, always_ranks = 2L, always_features = character()) {
  loadNamespace("mia")

  res <- mia::meltSE(se, add.row = TRUE, add.col = TRUE)
  all_ranks <- taxonomy_ranks(se)

  rank_names_show_always <- head(all_ranks, always_ranks)
  lowest_rank <- dplyr::last(all_ranks)

  # input assertion: exactly one value per sample and lowest rank
  stopifnot(identical(
    res |>
      dplyr::count(.data[[lowest_rank]], SampleID) |>
      dplyr::filter(n != 1L) |>
      nrow(),
    0L
  ))

  # zero count measurements slow things down and are often not required
  if (remove_zeros) {
    res <- dplyr::filter(res, counts > 0L)
  }

  # relative abundance
  res <-
    res |>
    group_by(SampleID) |>
    mutate(fraction = counts / sum(counts))

  # the "smart" part:
  res <- add_column(res, Feature = NA_character_)
  rank_names <- all_ranks
  for (rank in rev(rank_names)) {
    res <- tidyr::unite(res, ".lineage", all_of(rank_names), sep = " ‣ ", remove = FALSE)

    # find which lineages at current rank to still keep:
    if (rank %in% rank_names_show_always) {
      keep_linages <-
        res |>
        pull(.lineage) |>
        unique() |>
        sort()
    } else {
      keep_linages <-
        c(
          res |>
            dplyr::filter(is.na(Feature)) |>
            group_by(.lineage, SampleID) |>
            summarise(fraction = sum(fraction), .groups = "drop_last") |>
            dplyr::filter(fraction >= min_abundance) |>
            dplyr::count() |>
            dplyr::filter(n >= min_prevalence) |>
            pull(.lineage),
          res |>
            dplyr::filter(.data[[rank]] %in% always_features) |>
            pull(.lineage)
        ) |>
        unique() |>
        sort()
    }

    # set Feature column for all lineages to keep (which will be later used to agglomerate the data)
    # if the Feature column was already set, this means a previous rank already selected it and has priority
    res <- mutate(res, Feature = if_else(is.na(Feature) & .lineage %in% keep_linages, .lineage, Feature))

    # remove values from taxonomy column when feature is not kept:
    res <- mutate(res, "{rank}" := if_else(is.na(Feature), NA_character_, .data[[rank]]))

    # remove current rank: needed for unite() to work
    rank_names <- head(rank_names, -1L)
  }
  # clean up temporary column
  res <- res |> mutate(.lineage = NULL)

  # prepare feature data (average sequence length per feature)
  per_feature_data <-
    res |>
    group_by(Feature) |>
    summarise(
      sequence_length = sequence_length |> mean() |> signif(),
      decontam_p_value = decontam_p_value |> mean() |> signif()
    )

  # prepare sample data (yield / total read pairs per sample)
  per_sample_data <-
    res |>
    dplyr::filter(counts > 0L) |>
    group_by(SampleID) |>
    summarise(
      sample_count = sum(counts),
      n_features = n()
    )

  # agglomerate based on selected lineages
  res <-
    res |>
    group_by(across(all_of(c(
      all_ranks,
      "Feature",
      "SampleID",
      names(SummarizedExperiment::colData(se))
    )))) |>
    summarise(counts = sum(counts), fraction = sum(fraction), .groups = "drop")

  # merge all data ans relocate columns
  res <-
    res |>
    dplyr::left_join(per_feature_data, by = "Feature") |>
    dplyr::left_join(per_sample_data, by = "SampleID") |>
    relocate(Feature, SampleID, counts, fraction, all_of(all_ranks), all_of(names(per_feature_data)), all_of(names(per_sample_data)))

  res <- update_provenance(res, se, list(analysis = "bubble plot"))

  # save settings so it can be shown in plot caption
  attr(res, "min_abundance") <- min_abundance
  attr(res, "min_prevalence") <- min_prevalence
  attr(res, "remove_zeros") <- remove_zeros
  attr(res, "always_ranks") <- rank_names_show_always
  attr(res, "always_features") <- always_features

  # output assertion: exactly one value per sample and feature
  stopifnot(identical(
    res |>
      dplyr::count(Feature, SampleID) |>
      dplyr::filter(n != 1L) |>
      nrow(),
    0L
  ))

  # output assertion: the sum of all fractions per sample should be 1
  stopifnot(identical(
    res |>
      group_by(SampleID) |>
      summarise(fraction = sum(fraction)) |>
      dplyr::filter(round(fraction, 10L) != 1.0) |>
      nrow(),
    0L
  ))

  res
}

smart_bubble_plot <- function(
  orig_data,
  sample_label_from,
  facet_cols,
  theme,
  max_size = 10L
) {
  plot_data <-
    orig_data |>
    arrange(.data[[sample_label_from]], SampleID) |>
    mutate(
      Sample = str_c(
        "(",
        sample_count |> format(big.mark = " ", trim = TRUE),
        ") ",
        str_pad(SampleID, width = SampleID |> as.character() |> nchar() |> max(), side = "right", pad = " "),
        " ",
        str_pad(.data[[sample_label_from]], width = .data[[sample_label_from]] |> nchar() |> max(), pad = " ")
      ) |>
        fct_inorder(),
      Feature = Feature |>
        multi_taxa_trim() |>
        str_c(
          str_replace_na(str_c(" [", round(sequence_length, 1L), " bp]"), ""),
          if_else(!is.na(decontam_p_value), str_c(" {d=", sprintf("%.2f", decontam_p_value), "}"), "")
        ),
      fraction = fraction * 100.0 # RA (%)
    )

  caption <- str_c(
    "RA (relative abundance) shown for higher taxonomic ranks are exclusive of the RA for separately shown lower taxonomic ranks.",
    "\n",
    "Taxa are NOT collapsed if RA ≥ ",
    attr(plot_data, "min_abundance") * 100.0,
    "% in at least ",
    attr(plot_data, "min_prevalence"),
    " sample(s).",
    " ",
    "The highest rank(s) are never collapsed.",
    "\n",
    "Numbers in parentheses are total read pair counts (depth) per sample/library.",
    " ",
    "Bar chart is depth (filled bars) and number of features (open bars).",
    "\n",
    "Basepair numbers in brackets are (mean) feature nucleic acid sequence lengths.",
    "\n",
    "d-values are (mean) decontam results. Lower mean more likely a contamint."
  )

  main_plot <-
    ggplot(
      data = plot_data,
      mapping = aes(
        x = Sample,
        y = fct_rev(Feature),
        size = fraction,
        fill = Phylum,
        alpha = (!is.na(ASV_ID)) |> as.integer()
      )
    ) +
    geom_point(
      shape = 21L,
      colour = "black"
    ) +
    geom_text(
      data = plot_data |> dplyr::filter(fraction >= 1.0),
      mapping = aes(label = round(fraction)),
      size = 2.0,
      hjust = 0.35,
      vjust = 0.6
    ) +
    scale_size_area(
      max_size = max_size,
      breaks = 10L^(-9L:2L),
      labels = formatC
    ) +
    scale_alpha(
      range = c(0.6, 1.0),
      guide = "none"
    ) +
    facet_grid(
      cols = facet_cols,
      scales = "free",
      space = "free",
      labeller = label_wrap_gen(width = 10L)
    ) +
    labs(
      x = NULL,
      y = NULL,
      size = "RA (%)"
    ) +
    guides(fill = guide_legend(ncol = 1L)) +
    theme +
    theme(
      axis.text.x = element_text(angle = 90.0, hjust = 1.0, vjust = 0.5, family = "monospace")
    )

  main_plot <-
    main_plot |>
    update_provenance(orig_data) |>
    plot_titles()

  yield_data <-
    plot_data |>
    distinct(Sample, .keep_all = TRUE) |>
    pivot_longer(c(sample_count, n_features))

  scale_asv_count_by <- max(plot_data$sample_count) / max(plot_data$n_features)
  yield_data <- yield_data |> mutate(value = if_else(name == "n_features", value * scale_asv_count_by, value))

  yield_plot <-
    ggplot(
      data = yield_data,
      mapping = aes(
        x = Sample,
        y = value,
        fill = name |> fct_rev()
      )
    ) +
    geom_col(
      position = position_dodge2(padding = 0.25),
      colour = "darkgrey"
    ) +
    scale_y_reverse() +
    scale_fill_manual(
      values = c(sample_count = "darkgrey", n_features = "white"),
      guide = "none"
    ) +
    expand_limits(
      y = 0.0
    ) +
    facet_grid(
      cols = facet_cols,
      scales = "free",
      space = "free"
    ) +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = element_blank(),
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  if (yield_data |> has_name("Input_read_pairs")) {
    yield_plot <-
      yield_plot +
      geom_errorbar(
        mapping = aes(
          ymin = Input_read_pairs,
          ymax = Input_read_pairs,
          colour = name |> fct_rev()
        ),
        position = position_dodge2(padding = 0.25)
      ) +
      scale_color_manual(
        values = c(sample_count = "grey", n_features = "transparent"),
        guide = "none"
      )
  }

  p <-
    main_plot /
    yield_plot +
    plot_layout(
      heights = c(42.0, 1.0)
    ) +
    plot_annotation(
      caption = caption
    )

  p <- update_provenance(p, orig_data)

  n_samples <- plot_data |>
    dplyr::pull(Sample) |>
    vec_unique_count()
  n_features <- plot_data |>
    dplyr::pull(Feature) |>
    vec_unique_count()
  feature_names_width <- plot_data |>
    dplyr::pull(Feature) |>
    nchar() |>
    max()
  attr(p, "output") <- list(height = 8.0 + 0.15 * n_features, width = 0.05 * feature_names_width + 0.5 * n_samples)

  p
}
