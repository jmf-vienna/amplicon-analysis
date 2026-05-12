fill_unclassified <- function(
  data,
  ranks = any_of(c("Domain", "Kingdom")):ASV_ID,
  value = "unclassified",
  species_value = "sp.",
  parentheses = TRUE,
  verbose = TRUE
) {
  rank_names <- data |>
    dplyr::select({{ ranks }}) |>
    names()

  if (parentheses) {
    prefix <- "("
    suffix <- ")"
  } else {
    prefix <- suffix <- ""
  }

  prev_rank <- NA
  for (rank in rank_names) {
    if (is.na(prev_rank)) {
      data <- data |> mutate("{rank}" := .data[[rank]] |> replace_na(str_c(prefix, value, suffix)))
    } else {
      data <- data |>
        mutate(
          "{rank}" := coalesce(
            .data[[rank]],
            if_else(
              .data[[prev_rank]] |> str_detect(value),
              .data[[prev_rank]],
              .data[[prev_rank]] |> str_c(prefix, value, " ", prev_rank = _, suffix)
            )
          )
        )
    }

    if (!isFALSE(species_value) && prev_rank == "Genus" && rank == "Species") {
      data <- data |>
        mutate(
          Species = if_else(
            # if: Genus IS NOT unclassified & Species IS unclassified
            !str_detect(Genus, value) & str_detect(Species, value),
            # then: Genus sp.
            Genus |> str_c(prefix, prev_rank = _, " ", species_value, suffix),
            # else: leave as is
            Species
          )
        )
    }

    prev_rank <- rank
  }

  data
}

smart_agglomerate <- function(
  data,
  min_abundance = 1e-2,
  min_prevalence = 1,
  remove_zeros = TRUE,
  ranks = any_of(c("Domain", "Kingdom")):ASV_ID,
  always_ranks = 2,
  verbose = TRUE,
  always_features = c()
) {
  rank_names <- data |>
    dplyr::select({{ ranks }}) |>
    names()

  rank_names_show_always <- rank_names |> head(always_ranks)

  # Input assertion: exactly one value per sample and ASV/lowest rank
  stopifnot(identical(
    data |>
      dplyr::count(.data[[rank_names |> tail(1)]], Sample) |>
      dplyr::filter(n > 1) |>
      nrow(),
    0L
  ))

  # zero count measurements slow things down and are not always required
  if (remove_zeros) {
    data <- data |> dplyr::filter(Count > 0)
  }

  # relative abundance
  data <- data |>
    group_by(Sample) |>
    mutate(Fraction = Count / sum(Count))

  # set NAs to unclassified, in case this was not dealt with beforehand
  data <- data |> mutate(across(all_of(rank_names), \(x) replace_na(x, "(unclassified)")))

  # the "smart" part:
  data <- data |> add_column(Feature = NA_character_)
  for (rank in rev(rank_names)) {
    data <- data |> tidyr::unite(".lineage", all_of(rank_names), sep = " > ", remove = FALSE)

    # find which lineages at current rank to still keep:
    if (rank %in% rank_names_show_always) {
      keep_linages <-
        data |>
        pull(.lineage) |>
        unique() |>
        sort()
    } else {
      keep_linages <-
        c(
          data |>
            dplyr::filter(is.na(Feature)) |>
            group_by(.lineage, Sample) |>
            summarise(Fraction = sum(Fraction), .groups = "drop_last") |>
            dplyr::filter(Fraction >= min_abundance) |>
            dplyr::count() |>
            dplyr::filter(n >= min_prevalence) |>
            pull(.lineage),
          data |>
            dplyr::filter(.data[[rank]] %in% always_features) |>
            pull(.lineage)
        ) |>
        sort()
    }
    # keep_linages |> cat(sep = "\n")

    # set Feature column for all lineages to keep (which will be later used to agglomerate the data)
    # if the Feature column was already set, this means a previous rank already selected it and has priority
    data <- data |> mutate(Feature = if_else(is.na(Feature) & .lineage %in% keep_linages, .lineage, Feature))

    # remove values from taxonomy column when feature is not kept:
    data <- data |> mutate("{rank}" := if_else(is.na(Feature), NA_character_, .data[[rank]]))

    # remove current rank - needed for unite() to work
    rank_names <- head(rank_names, -1)
  }
  # clean up temporary column
  data <- data |> mutate(.lineage = NULL)

  # prepare feature data (average sequence length per feature)
  per_feature_data <-
    data |>
    group_by(Feature) |>
    summarise(Sequence_length = Sequence_length |> mean())

  # prepare sample data (yield / total read pairs per sample)
  per_sample_data <-
    data |>
    dplyr::filter(Count > 0) |>
    group_by(Sample) |>
    summarise(
      Sample_count = sum(Count),
      ASVs = n()
    )

  # agglomerate based on selected lineages
  data <-
    data |>
    group_by(across(
      c({{ ranks }}, Sample) |
        !any_of(c(
          "Count",
          "Fraction",
          "Sequence",
          "Sequence_length",
          "sequence",
          "sequence_length",
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
    summarise(Count = sum(Count), Fraction = sum(Fraction), .groups = "drop")

  # merge all data ans relocate columns
  data <-
    data |>
    dplyr::left_join(per_feature_data, by = "Feature") |>
    dplyr::left_join(per_sample_data, by = "Sample") |>
    relocate(Feature, Sample, Group, Count, Fraction, {{ ranks }}, Sequence_length)

  # save some settings so it can be shown in plot caption
  attr(data, "min_abundance") <- min_abundance
  attr(data, "min_prevalence") <- min_prevalence

  # Output assertion: exactly one value per sample and feature
  stopifnot(identical(
    data |>
      dplyr::count(Feature, Sample) |>
      dplyr::filter(n > 1) |>
      nrow(),
    0L
  ))

  # Output assertion: the sum of all fractions per sample should be 1
  stopifnot(identical(
    data |>
      group_by(Sample) |>
      summarise(Fraction = sum(Fraction)) |>
      dplyr::filter(round(Fraction, 10) != 1) |>
      nrow(),
    0L
  ))

  data
}

smart_agglomerate_bubble_plot <- function(
  data,
  title = waiver(),
  subtitle = waiver(),
  caption = NA,
  colour = Phylum,
  facets = vars(Group),
  facet_labeller = label_wrap_gen(width = 10),
  max_size = 10,
  add_sequence_length = TRUE,
  trim_multi_taxa = FALSE,
  verbose = TRUE
) {
  orig_data <- data

  require(patchwork)

  data <- data |>
    arrange(Sample) |>
    mutate(
      # pad sample name, and move the spaces between ID and user sample name (or wherever the first space is) + add sample count
      Sample = Sample |>
        str_pad(data |> pull(Sample) |> nchar() |> max(), pad = " ") |>
        str_replace("^( +)([^ ]+)", "\\2\\1") |>
        # Unicode 2009 is https://en.wikipedia.org/wiki/Thin_space
        # str_c("{", ASVs |> format(big.mark = "\u2009"), "} ", sample = _) |>
        str_c("(", Sample_count |> format(big.mark = " ", trim = TRUE), ") ", sample = _) |>
        fct_inorder(),
      # RA (%)
      Fraction = Fraction * 100
    )

  if (add_sequence_length) {
    data <- data |>
      mutate(
        Feature = str_c(Feature, str_replace_na(str_c(" [", Sequence_length |> round(1), " bp]"), ""))
      )
  }

  if (data |> has_name("Feature_notes")) {
    data <- data |>
      mutate(
        Feature = str_c(Feature, str_replace_na(str_c(" {", Feature_notes, "}"), ""))
      )
  }

  if (trim_multi_taxa) {
    data <- data |>
      mutate(
        # for multi species trimming: "Genus many/different/possible/species/names" -> "Genus many/…/names"
        Feature = Feature |> str_replace_all("/[a-z\\./]+/", "/…/")
      )
  }

  caption <- str_c(
    "RA (relative abundance) shown for higher taxonomic ranks are exclusive of the RA for separately shown lower taxonomic ranks.",
    "\n",
    "Taxa are NOT collapsed if RA ≥ ",
    attr(data, "min_abundance") * 100,
    "% in at least ",
    attr(data, "min_prevalence"),
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

  if (verbose) {
    message("Caption:")
    cat(caption)

    message("Features:")
    data |>
      pull(Feature) |>
      unique() |>
      sort() |>
      cat(sep = "\n")
  }

  main_plot <- ggplot(
    data = data,
    mapping = aes(
      x = Sample,
      y = Feature |> fct_rev(),
      size = Fraction,
      fill = {{ colour }},
      alpha = (!is.na(ASV_ID)) |> as.integer()
    )
  ) +
    geom_point(
      shape = 21,
      colour = "black"
    ) +
    geom_text(
      data = data |> dplyr::filter(Fraction >= 1),
      mapping = aes(label = round(Fraction)),
      size = 2,
      hjust = 0.35,
      vjust = 0.6,
    ) +
    scale_size_area(
      max_size = max_size,
      breaks = 10^(-9:2),
      labels = formatC
    ) +
    scale_alpha(
      range = c(0.6, 1),
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
    guides(fill = guide_legend(ncol = 1)) +
    theme(
      axis.text = ggplot2::element_text(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.subtitle = element_text(colour = "blue"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = "monospace")
    )

  yield_data <-
    data |>
    distinct(Sample, .keep_all = TRUE) |>
    pivot_longer(c(Sample_count, ASVs))

  scale_ASV_count_by <- max(data$Sample_count) / max(data$ASVs)
  yield_data <- yield_data |> mutate(value = if_else(name == "ASVs", value * scale_ASV_count_by, value))

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
      y = 0
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
      heights = c(42, 1)
    ) +
    plot_annotation(
      caption = caption,
    )

  p <- update_provenance(p, orig_data)

  n_samples <- data |> dplyr::pull(Sample) |> vec_unique_count()
  n_features <- data |> dplyr::pull(Feature) |> vec_unique_count()
  feature_names_width <- data |> dplyr::pull(Feature) |> nchar() |> max()
  attr(p, "output") <- list(height = 8.0 + 0.15 * n_features, width = 0.05 * feature_names_width + 0.5 * n_samples)

  p
}
