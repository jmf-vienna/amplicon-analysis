library(targets)

jmf::quiet()
tar_option_set(format = "qs")

tar_config_get("script") |>
  fs::path_dir() |>
  fs::path("R") |>
  fs::dir_ls() |>
  purrr::walk(source)

list(
  tar_target(pipeline_version, "0.1.0"),
  tar_target(theme, ggplot_theme()),

  # config:
  tar_target(results_dir_name, fs::dir_ls("Results", type = "dir") |> head(1L)),
  tar_target(project_name, results_dir_name |> stringr::str_extract(jmf::jmf_project_id_regex())),
  tar_target(plots_dir_name, "plots"),

  # assay data (counts):
  tar_target(counts_file_name, fs::path(results_dir_name, "DADA2_counts.tsv")),
  tar_target(counts_file, counts_file_name, format = "file"),
  tar_target(counts, readr::read_tsv(counts_file)),
  tar_target(assay_data, make_assay_data(counts)),
  tar_target(debug_assay_data, str(assay_data)),

  # column data (libraries, samples):
  tar_target(libraries_file_name, fs::path("Metadata", "Libraries.tsv")),
  tar_target(libraries_file, libraries_file_name, format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file)),
  tar_target(libraries_col_data, make_col_data(libraries)),
  tar_target(debug_col_data, print(libraries_col_data)),

  # row data (taxonomy):
  tar_target(taxonomy_file_name, fs::path(results_dir_name, "DADA2_ASVs.rRNA_SSU.SILVA_reference.DADA2_classified.tsv")),
  tar_target(taxonomy_file, taxonomy_file_name, format = "file"),
  tar_target(taxonomy, readr::read_tsv(taxonomy_file)),
  tar_target(row_data, make_row_data(taxonomy)),
  tar_target(debug_row_data, print(row_data)),

  # SummarizedExperiment, raw:
  tar_target(se_raw, se(assay_data, libraries_col_data, row_data, list(project = project_name, state = "raw"))),
  tar_target(debug_se_raw, print(se_raw)),

  # ordination:
  tar_target(ps, as_phyloseq(se_raw)),
  tar_target(ps_distance, calulcate_distance(ps)),
  tar_target(ps_ordination, calulcate_ordination(ps_distance)),
  tar_target(ordination_plot, plot_ordination(ps_ordination, "Sequencing_date", theme = theme)),
  tar_target(ordination_plot_file, save_plot(ordination_plot, plots_dir_name), format = "file")
)
