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
  glob <- glue::glue("*.{reference}_reference.{classifier}_classified.tsv")

  file <- fs::dir_ls(path, glob = glob)

  if (vec_is_empty(file)) {
    cli::cli_alert_danger("failed finding taxonomy file using {.val {glob}}")
    return(invisible())
  }

  if (length(file) > 1L) {
    file <- file |> head(1L)
  }

  cli::cli_alert("found taxonomy table at {.file {file}}")
  file
}

find_features_info_file <- function(path) {
  file <- fs::dir_ls(path, glob = "*ASVs.tsv")
  cli::cli_alert("found features info table at {.file {file}}")
  file
}
