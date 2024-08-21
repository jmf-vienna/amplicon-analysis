ggplot_theme <- function() {
  ggplot2::theme_gray() + ggplot2::theme(
    plot.background = ggplot2::element_blank(),
    text = ggplot2::element_text(family = "Noto Sans"),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.ticks = ggplot2::element_line(colour = "black")
  )
}

save_plot <- function(plot, file_name) {
  svglite::svglite(
    fs::path(file_name, ext = "svg"),
    web_fonts = list("https://fonts.googleapis.com/css?family=Noto%20Sans")
  )
  print(plot)
  dev.off()
}
