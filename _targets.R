library(targets)

jmf::quiet()
tar_option_set(format = "qs")

tar_config_get("script") |>
  fs::path_dir() |>
  fs::path("R") |>
  fs::dir_ls() |>
  purrr::walk(source)

list(
  # Column data:
  tar_target(libraries_file_name, fs::path("Metadata", "Libraries.tsv")),
  tar_target(libraries_file, libraries_file_name, format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file)),
  tar_target(col_data, make_col_data(libraries)),
  tar_target(debug_col_data, print(col_data))
)
