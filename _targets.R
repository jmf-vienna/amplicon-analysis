library(targets)

jmf::quiet()
tar_option_set(format = "qs")

tar_config_get("script") |>
  fs::path_dir() |>
  fs::path("R") |>
  fs::dir_ls() |>
  purrr::walk(source)

list(
  tar_target(results_dir_name, fs::dir_ls("Results", type = "dir") |> head(1)),

  # Count assay data:
  tar_target(counts_file_name, fs::path(results_dir_name, "DADA2_counts.tsv")),
  tar_target(counts_file, counts_file_name, format = "file"),
  tar_target(counts, readr::read_tsv(counts_file)),
  tar_target(assay_data, make_assay_data(counts)),
  tar_target(debug_assay_data, print(assay_data)),

  # Column data:
  tar_target(libraries_file_name, fs::path("Metadata", "Libraries.tsv")),
  tar_target(libraries_file, libraries_file_name, format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file)),
  tar_target(col_data, make_col_data(libraries)),
  tar_target(debug_col_data, print(col_data)),

  # Row data:
  tar_target(taxonomy_file_name, fs::path(results_dir_name, "DADA2_ASVs.rRNA_SSU.SILVA_reference.DADA2_classified.tsv")),
  tar_target(taxonomy_file, taxonomy_file_name, format = "file"),
  tar_target(taxonomy, readr::read_tsv(taxonomy_file)),
  tar_target(row_data, make_row_data(taxonomy)),
  tar_target(debug_row_data, print(row_data))
)
