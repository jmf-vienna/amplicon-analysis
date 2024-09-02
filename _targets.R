library(targets)

jmf::quiet()
options(warn = 2L)
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

  # config > paths:
  tar_target(data_dir_name, config[["path"]][["data"]]),
  tar_target(results_dir_name, config[["path"]][["results"]]),
  tar_target(plots_dir_name, config[["path"]][["plots"]]),

  # config > refinement:
  tar_target(desirables, config[["filter"]][["desirable"]]),
  tar_target(undesirables, config[["filter"]][["undesirable"]]),
  tar_target(yield_min, config[["filter"]][["yield"]][["min"]]),

  # config > sample data column names:
  tar_target(sample_label_from, config[["annotation"]][["sample"]]),
  tar_target(variable_of_interest, config[["analyse"]][["category"]]),

  # column data > samples:
  tar_target(samples_file, "Samples.tsv", format = "file"),
  tar_target(samples, readr::read_tsv(samples_file)),
  tar_target(samples_col_data, make_col_data(list(samples))),

  # column data > libraries:
  tar_target(libraries_file, "Libraries.tsv", format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file)),
  tar_target(libraries_col_data, make_col_data(list(libraries, samples))),

  # assay data (counts):
  tar_target(counts_file, find_counts_file(data_dir_name), format = "file"),
  tar_target(counts, readr::read_tsv(counts_file)),
  tar_target(assay_data, make_assay_data(counts)),

  # row data (taxonomy):
  tar_target(taxonomy_file, find_taxonomy_file(data_dir_name, config[["taxonomy"]]), format = "file"),
  tar_target(taxonomy, readr::read_tsv(taxonomy_file)),
  tar_target(row_data, make_row_data(taxonomy)),

  # SummarizedExperiment:
  tar_target(base_provenance, list(project = config[["project"]][["name"]], gene = config[["gene"]][["name"]])),

  # SummarizedExperiment > libraries > raw:
  tar_target(se_libs_raw_provenance, modifyList(base_provenance, list(stage = "libraries", state = "raw"))),
  tar_target(se_libs_raw, make_se(assay_data, libraries_col_data, row_data, se_libs_raw_provenance)),

  # SummarizedExperiment > samples > raw:
  tar_target(se_raw, merge_cols(se_libs_raw, samples |> names() |> head(1L), samples |> names(), list(stage = "samples"))),

  # SummarizedExperiment > samples > refined:
  tar_target(
    se_refined,
    se_raw |>
      keep_desirable_features(desirables) |>
      filter_undesirable_features(undesirables) |>
      update_provenance(se_raw, list(state = "refined")) |>
      tidy()
  ),
  tar_target(
    se_deep,
    se_refined |>
      filter_samples_by_sum(yield_min) |>
      tidy()
  ),

  # SummarizedExperiment > all:
  tar_target(se, list(se_libs_raw, se_raw, se_refined, se_deep), iteration = "list"),
  tar_target(se_flat_file, export_flattened(se, results_dir_name), format = "file", pattern = map(se)),
  tar_target(ps, as_phyloseq(se), pattern = map(se)),

  # ordination:
  tar_target(ps_distance, calulcate_distance(ps), pattern = map(ps)),
  tar_target(ps_ordination, calulcate_ordination(ps_distance), pattern = map(ps_distance)),
  tar_target(
    ordination_plot,
    plot_ordination(ps_ordination, variable_of_interest, sample_label_from, theme),
    pattern = cross(ps_ordination, variable_of_interest)
  ),
  tar_target(ordination_plot_file, save_plot(ordination_plot, plots_dir_name), format = "file", pattern = map(ordination_plot))
)
