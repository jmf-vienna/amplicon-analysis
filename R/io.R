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

force_valid_file_name <- function(x) {
  x |> stringr::str_replace_all("[^a-zA-Z0-9-]", "_")
}

prepare_export <- function(file) {
  file |>
    fs::path_dir() |>
    fs::dir_create(mode = Sys.getenv("DIR_CREATE_MODE", "u=rwx,go=rx"))
}

write_rds <- function(x, file) {
  prepare_export(file)
  saveRDS(x, file)
  cli::cli_alert("RDS saved to {.file {file}}")
  file
}

write_qs2 <- function(x, file) {
  prepare_export(file)
  qs2::qs_save(x, file)
  cli::cli_alert("QS2 saved to {.file {file}}")
  file
}

write_tsv <- function(x, file, na = "invisible") {
  prepare_export(file)

  rlang::arg_match0(na, c("invisible", "explicit"))
  na <- dplyr::case_when(
    na == "invisible" ~ "",
    na == "explicit" ~ "NA"
  )

  readr::write_tsv(x, file, na = na)
  cli::cli_alert("table saved to {.file {file}}")
  file
}

write_text <- function(x, file) {
  prepare_export(file)
  cat(x, file = file)
  cli::cli_alert("text saved to {.file {file}}")
  file
}
