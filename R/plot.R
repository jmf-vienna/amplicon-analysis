save_plot <- function(plot, file_name) {
  svglite::svglite(fs::path(file_name, ext = "svg"))
  print(plot)
  dev.off()
}
