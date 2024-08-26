ggplot_theme <- function() {
  ggplot2::theme_gray() + ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    text = ggplot2::element_text(family = "Noto Sans"),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.ticks = ggplot2::element_line(colour = "black")
  )
}

save_plot <- function(plot, dir_name) {
  file <- fs::path(dir_name, plot |> as_file_name(), ext = "svg")

  prepare_export(file)

  svglite::svglite(
    file,
    web_fonts = list("https://fonts.googleapis.com/css?family=Noto%20Sans")
  )
  print(plot)
  dev.off()

  invisible(file)
}
