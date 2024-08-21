prepare_export <- function(file) {
  file |>
    fs::path_dir() |>
    fs::dir_create()
}
