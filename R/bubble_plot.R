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

smart_agglomerate <- function(se, min_abundance = 0.01, min_prevalence = 2L, remove_zeros = TRUE, always_ranks = 2L, always_features = NULL) {
  loadNamespace("mia")

  res <- mia::meltSE(se, add.row = TRUE, add.col = TRUE)

  rank_names <- all_ranks <- taxonomy_ranks(se)

  rank_names_show_always <- all_ranks |> head(always_ranks)

  # Input assertion: exactly one value per sample and ASV/lowest rank
  stopifnot(identical(
    res |>
      dplyr::count(.data[[tail(all_ranks, 1L)]], SampleID) |>
      dplyr::filter(n > 1L) |>
      nrow(),
    0L
  ))

  # zero count measurements slow things down and are not always required
  if (remove_zeros) {
    res <- res |> dplyr::filter(counts > 0L)
  }

  # relative abundance
  res <-
    res |>
    group_by(SampleID) |>
    mutate(Fraction = counts / sum(counts))

  # the "smart" part:
  res <- res |> add_column(Feature = NA_character_)
  for (rank in rev(rank_names)) {
    res <- res |> tidyr::unite(".lineage", all_of(rank_names), sep = " > ", remove = FALSE)

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
            summarise(Fraction = sum(Fraction), .groups = "drop_last") |>
            dplyr::filter(Fraction >= min_abundance) |>
            dplyr::count() |>
            dplyr::filter(n >= min_prevalence) |>
            pull(.lineage),
          res |>
            dplyr::filter(.data[[rank]] %in% always_features) |>
            pull(.lineage)
        ) |>
        sort()
    }

    # set Feature column for all lineages to keep (which will be later used to agglomerate the data)
    # if the Feature column was already set, this means a previous rank already selected it and has priority
    res <- res |> mutate(Feature = if_else(is.na(Feature) & .lineage %in% keep_linages, .lineage, Feature))

    # remove values from taxonomy column when feature is not kept:
    res <- res |> mutate("{rank}" := if_else(is.na(Feature), NA_character_, .data[[rank]]))

    # remove current rank - needed for unite() to work
    rank_names <- head(rank_names, -1L)
  }
  # clean up temporary column
  res <- res |> mutate(.lineage = NULL)

  # prepare feature data (average sequence length per feature)
  per_feature_data <-
    res |>
    group_by(Feature) |>
    summarise(sequence_length = sequence_length |> mean())

  # prepare sample data (yield / total read pairs per sample)
  per_sample_data <-
    res |>
    dplyr::filter(counts > 0L) |>
    group_by(SampleID) |>
    summarise(
      Sample_count = sum(counts),
      ASVs = n()
    )

  # agglomerate based on selected lineages
  res <-
    res |>
    group_by(across(
      c(all_of(all_ranks), SampleID) |
        !any_of(c(
          "counts",
          "Fraction",
          "Sequence",
          "sequence_length",
          "FeatureID",
          "sequence",
          "decontam_p_value",
          "Orientation",
          "quality_min_eepm",
          "sha1",
          "sha1base36",
          "Lineage",
          "bootstrap_Kingdom",
          "bootstrap_Domain",
          "bootstrap_Phylum",
          "bootstrap_Class",
          "bootstrap_Order",
          "bootstrap_Family",
          "bootstrap_Genus",
          "bootstrap_Species",
          "BLAST_percent_identity",
          "Strain"
        ))
    )) |>
    summarise(counts = sum(counts), Fraction = sum(Fraction), .groups = "drop")

  # merge all data ans relocate columns
  res <-
    res |>
    dplyr::left_join(per_feature_data, by = "Feature") |>
    dplyr::left_join(per_sample_data, by = "SampleID") |>
    relocate(Feature, SampleID, Group, counts, Fraction, all_of(all_ranks), sequence_length)

  # save some settings so it can be shown in plot caption
  attr(res, "min_abundance") <- min_abundance
  attr(res, "min_prevalence") <- min_prevalence

  # Output assertion: exactly one value per sample and feature
  stopifnot(identical(
    res |>
      dplyr::count(Feature, SampleID) |>
      dplyr::filter(n > 1L) |>
      nrow(),
    0L
  ))

  # Output assertion: the sum of all fractions per sample should be 1
  stopifnot(identical(
    res |>
      group_by(SampleID) |>
      summarise(Fraction = sum(Fraction)) |>
      dplyr::filter(round(Fraction, 10L) != 1L) |>
      nrow(),
    0L
  ))

  res |> update_provenance(se, list(analysis = "bubble plot"))
}

smart_bubble_plot <- function(
  data,
  sample_label_from,
  title = waiver(),
  subtitle = waiver(),
  caption = NA,
  colour = Phylum,
  facets = vars(Group),
  facet_labeller = label_wrap_gen(width = 10L),
  max_size = 10L,
  add_sequence_length = TRUE,
  trim_multi_taxa = FALSE
) {
  orig_data <- data

  plot_data <-
    data |>
    mutate(
      Sample = str_c(SampleID, " ", .data[[sample_label_from]])
    )

  plot_data <-
    plot_data |>
    arrange(Sample) |>
    mutate(
      # pad sample name, and move the spaces between ID and user sample name (or wherever the first space is) + add sample count
      Sample = Sample |>
        str_pad(plot_data |> pull(Sample) |> nchar() |> max(), pad = " ") |>
        str_replace("^( +)([^ ]+)", "\\2\\1") |>
        # Unicode 2009 is https://en.wikipedia.org/wiki/Thin_space
        # str_c("{", ASVs |> format(big.mark = "\u2009"), "} ", sample = _) |>
        str_c("(", Sample_count |> format(big.mark = " ", trim = TRUE), ") ", sample = _) |>
        fct_inorder(),
      # RA (%)
      Fraction = Fraction * 100.0
    )

  if (add_sequence_length) {
    plot_data <- plot_data |>
      mutate(
        Feature = str_c(Feature, str_replace_na(str_c(" [", sequence_length |> round(1L), " bp]"), ""))
      )
  }

  if (plot_data |> has_name("Feature_notes")) {
    plot_data <- plot_data |>
      mutate(
        Feature = str_c(Feature, str_replace_na(str_c(" {", Feature_notes, "}"), ""))
      )
  }

  if (trim_multi_taxa) {
    plot_data <- plot_data |>
      mutate(
        # for multi species trimming: "Genus many/different/possible/species/names" -> "Genus many/…/names"
        Feature = Feature |> str_replace_all("/[a-z\\./]+/", "/…/")
      )
  }

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
    "Bar chart is depth (filled bars) and number of ASVs (open bars).",
    if_else(add_sequence_length, str_c("\n", "Basepair numbers in brackets are (mean) ASV nucleic acid sequence lengths."), ""),
    str_replace_na(str_c("\n", caption), "")
  )

  main_plot <- ggplot(
    data = plot_data,
    mapping = aes(
      x = Sample,
      y = Feature |> fct_rev(),
      size = Fraction,
      fill = {{ colour }},
      alpha = (!is.na(ASV_ID)) |> as.integer()
    )
  ) +
    geom_point(
      shape = 21L,
      colour = "black"
    ) +
    geom_text(
      data = plot_data |> dplyr::filter(Fraction >= 1.0),
      mapping = aes(label = round(Fraction)),
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
      cols = facets,
      scales = "free",
      space = "free",
      labeller = facet_labeller
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      size = "RA (%)"
    ) +
    guides(fill = guide_legend(ncol = 1L)) +
    theme(
      axis.text = ggplot2::element_text(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.subtitle = element_text(colour = "blue"),
      axis.text.x = element_text(angle = 90.0, hjust = 1.0, vjust = 0.5, family = "monospace")
    )

  yield_data <-
    plot_data |>
    distinct(Sample, .keep_all = TRUE) |>
    pivot_longer(c(Sample_count, ASVs))

  scale_asv_count_by <- max(plot_data$Sample_count) / max(plot_data$ASVs)
  yield_data <- yield_data |> mutate(value = if_else(name == "ASVs", value * scale_asv_count_by, value))

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
      colour = "blue"
    ) +
    scale_y_reverse() +
    scale_fill_manual(
      values = c(Sample_count = "blue", ASVs = "white"),
      guide = "none"
    ) +
    expand_limits(
      y = 0.0
    ) +
    facet_grid(
      cols = facets,
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
        values = c(Sample_count = "grey", ASVs = "transparent"),
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
