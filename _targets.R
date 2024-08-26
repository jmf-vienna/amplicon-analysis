library(targets)

jmf::quiet()
tar_option_set(format = "qs")

tar_config_get("script") |>
  fs::path_dir() |>
  fs::path("R") |>
  fs::dir_ls(glob = "*.R") |>
  purrr::walk(source)

if (fs::dir_exists("R")) {
  fs::path("R") |>
    fs::dir_ls(glob = "*.R") |>
    purrr::walk(source)
}

list(
  tar_target(pipeline_version, "0.1.0"),
  tar_target(theme, ggplot_theme()),

  # config:
  tar_target(config_file, "config.yaml", format = "file"),
  tar_target(config, config::get(file = config_file)),
  tar_target(debug.config, str(config)),
  tar_target(data_dir_name, config[["path"]][["data"]]),
  tar_target(plots_dir_name, config[["path"]][["plots"]]),
  tar_target(categories, config[["analyse"]][["categories"]]),

  # column data > samples:
  tar_target(samples_file, "Samples.tsv", format = "file"),
  tar_target(samples, readr::read_tsv(samples_file)),
  tar_target(samples_col_data, make_col_data(list(samples))),
  tar_target(debug.samples_col_data, print(samples_col_data)),

  # column data > libraries:
  tar_target(libraries_file, "Libraries.tsv", format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file)),
  tar_target(libraries_col_data, make_col_data(list(libraries, samples))),
  tar_target(debug.libraries_col_data, print(libraries_col_data)),

  # assay data (counts):
  tar_target(counts_file, find_counts_file(data_dir_name), format = "file"),
  tar_target(counts, readr::read_tsv(counts_file)),
  tar_target(assay_data, make_assay_data(counts)),
  tar_target(debug.assay_data, str(assay_data)),

  # row data (taxonomy):
  tar_target(taxonomy_file, find_taxonomy_file(data_dir_name, config[["taxonomy"]]), format = "file"),
  tar_target(taxonomy, readr::read_tsv(taxonomy_file)),
  tar_target(row_data, make_row_data(taxonomy)),
  tar_target(debug.row_data, print(row_data)),

  # SummarizedExperiment, libraries, raw:
  tar_target(base_provenance, list(project = config[["project"]][["name"]], gene = config[["gene"]][["name"]])),
  tar_target(se_libs_raw_provenance, modifyList(base_provenance, list(stage = "libraries", state = "raw"))),
  tar_target(se_libs_raw, se(assay_data, libraries_col_data, row_data, se_libs_raw_provenance)),
  tar_target(debug.se_libs_raw, print(se_libs_raw)),

  # ordination:
  tar_target(ps, as_phyloseq(se_libs_raw)),
  tar_target(ps_distance, calulcate_distance(ps)),
  tar_target(ps_ordination, calulcate_ordination(ps_distance)),
  tar_target(ordination_plot, plot_ordination(ps_ordination, categories |> head(1L), theme = theme)),
  tar_target(ordination_plot_file, save_plot(ordination_plot, plots_dir_name), format = "file")
)
