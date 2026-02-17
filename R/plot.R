ggplot_theme <- function() {
  ggplot2::theme_gray() +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      text = ggplot2::element_text(family = font_family()),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

font_family <- function() {
  "Noto Sans"
}

save_plot <- function(x, dir_name, ...) {
  if (is.null(x)) {
    return(invisible())
  }

  file <- fs::path(dir_name, provenance_as_file_name(x), ext = "svg")
  prepare_export(file)

  svglite::svglite(
    file,
    web_fonts = list("https://fonts.googleapis.com/css?family=Noto%20Sans"),
    ...
  )
  print(x)
  dev.off()

  cli::cli_alert("plot saved to {.file {file}}")
  invisible(file)
}
