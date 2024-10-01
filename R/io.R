find_counts_file <- function(path) {
  file <- fs::dir_ls(path, glob = "*counts.tsv")
  cli::cli_alert("found counts table at {.file {file}}")
  file
}

find_libraries_summary_file <- function(path) {
  file <- fs::dir_ls(path, glob = "*libraries.tsv")
  cli::cli_alert("found libraries summary table at {.file {file}}")
  file
}

find_taxonomy_file <- function(path, params) {
  reference <- params |> purrr::pluck("reference", .default = "*")
  classifier <- params |> purrr::pluck("classifier", .default = "*")
  file <- fs::dir_ls(path, glob = glue::glue("*.{reference}_reference.{classifier}_classified.tsv")) |> head(1L)
  cli::cli_alert("found taxonomy table at {.file {file}}")
  file
}

find_features_info_file <- function(path) {
  file <- fs::dir_ls(path, glob = "*ASVs.tsv")
  cli::cli_alert("found features info table at {.file {file}}")
  file
}

force_valid_file_name <- function(x) {
  x |> stringr::str_replace_all("[^a-zA-Z0-9-]", "_")
}

prepare_export <- function(file) {
  file |>
    fs::path_dir() |>
    fs::dir_create(mode = Sys.getenv("DIR_CREATE_MODE", "u=rwx,go=rx"))
}

write_tsv <- function(x, file) {
  prepare_export(file)
  readr::write_tsv(x, file, na = "")
  cli::cli_alert("table saved to {.file {file}}")
  file
}

write_text <- function(x, file) {
  prepare_export(file)
  cat(x, file = file)
  cli::cli_alert("text saved to {.file {file}}")
  file
}
