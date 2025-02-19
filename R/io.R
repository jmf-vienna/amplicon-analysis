force_valid_file_name <- function(x) {
  x |> stringr::str_replace_all("[^a-zA-Z0-9-]", "_")
}

find_one_file <- function(path, glob, verbose = TRUE) {
  res <- dir_ls(path, type = "file", glob = glob)

  if (!is_string(res)) {
    cli_abort("Expected exactly one file matching the pattern {.arg {glob}} in {.path {path}}, but found {?none/}{.file {res}} instead!")
  }

  if (verbose) {
    cli_alert("found {.file {res}}")
  }

  res
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
