find_counts_file <- function(path) {
  file <- fs::dir_ls(path, glob = "*counts.tsv")
  cli::cli_alert("found counts table at {.file {file}}")
  file
}

find_taxonomy_file <- function(path, config) {
  file <- fs::dir_ls(path, glob = glue::glue_data(config, "*.{reference}_reference.{classfier}_classified.tsv"))
  cli::cli_alert("found taxonomy table at {.file {file}}")
  file
}

prepare_export <- function(file) {
  file |>
    fs::path_dir() |>
    fs::dir_create()
}
