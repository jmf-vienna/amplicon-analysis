calulcate_distance <- function(ps, distance = "aitchison") {
  ps |>
    microViz::dist_calc(distance) |>
    update_provenance(ps, list(distance = distance))
}

calulcate_ordination <- function(ps) {
  ps_new <-
    ps |>
    microViz::ord_calc()

  ps_new |>
    update_provenance(ps, list(
      ordination = microViz::info_get(ps_new)[["ord_info"]][["method"]]
    ))
}

plot_ordination <- function(ps, group, point_label, theme) {
  v <- ps |>
    phyloseq::sample_data() |>
    dplyr::pull({{ group }})
  vi <- variable_info(v)

  plot <-
    ps |>
    microViz::ord_plot(
      # NOTE: color and colour behave differently in microViz
      colour = group,
      fill = group
    ) + theme + ggplot2::theme(
      aspect.ratio = 1.0
    )

  if (vi[["gte3"]][["unique"]] > 1L) {
    plot <- plot + microViz::stat_chull(
      mapping = ggplot2::aes(colour = .data[[group]], fill = .data[[group]]),
      alpha = 0.1
    )
  }

  plot <- plot +
    ggrepel::geom_text_repel(
      mapping = ggplot2::aes(
        label = .data[[point_label]]
      ),
      max.overlaps = 100L,
      seed = 0L,
      size = 2.5,
      family = "Noto Sans"
    )

  plot |>
    plot_titles(ps) |>
    update_provenance(ps, list(
      aesthetics = list(color = group)
    ))
}
