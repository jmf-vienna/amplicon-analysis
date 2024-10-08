library(purrr)
library(targets)

jmf::quiet()
options(warn = 2L)
targets::tar_option_set(
  packages = c("cli", "stringr", "vctrs"),
  format = "qs",
  iteration = "list"
)

targets::tar_config_get("script") |>
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
  tar_target(pipeline_version, get_pipeline_version(), cue = targets::tar_cue(mode = "always")),
  tar_target(theme, ggplot_theme()),

  # config ----
  tar_target(config_file, Sys.getenv("R_CONFIG_FILE", "config.yaml"), format = "file"),
  tar_target(config, config::get(config = Sys.getenv("TAR_PROJECT", "default"), file = config_file)),
  tar_target(project_name, config |> pluck("project", "name", .default = "SOME PROJECT")),
  tar_target(gene_name, config |> pluck("gene", "name", .default = "SOME GENE")),
  tar_target(base_provenance, list(project = project_name, gene = gene_name)),

  ## config > io ----
  tar_target(data_dir_name, config |> pluck("path", "data", .default = "data")),
  tar_target(results_dir_name, config |> pluck("path", "results", .default = "results")),
  tar_target(plots_dir_name, config |> pluck("path", "plots", .default = "plots")),
  tar_target(file_prefix, base_provenance |> as_file_name()),

  ## config > refinement ----
  tar_target(negative_controls, config |> pluck("filter", "negative controls", .default = character())),
  tar_target(desirables, config |> pluck("filter", "desirable", .default = list())),
  tar_target(undesirables, config |> pluck("filter", "undesirable", .default = list())),
  tar_target(yield_min, config |> pluck("filter", "yield", "min", .default = 0L)),
  tar_target(yield_max, config |> pluck("filter", "yield", "max", .default = Inf)),
  tar_target(decontam_threshold, config |> pluck("filter", "decontam", .default = 0.0)),

  ## config > sample data column names ----
  tar_target(sample_label_from, config |> pluck("annotation", "sample", "variable name")),
  tar_target(variable_of_interest, config |> pluck("analyze", "category", .default = "Category")),

  ## config > limits ----
  tar_target(limits, list(
    sample = config |> pluck("annotation", "sample", "limit", .default = 10L),
    variable_of_interest = config |> pluck("annotation", "category", "limit", .default = 10L)
  )),

  # column data > samples ----
  tar_target(samples_file, "Samples.tsv", format = "file"),
  tar_target(samples, readr::read_tsv(samples_file)),
  tar_target(samples_col_data, make_col_data(list(samples))),

  # column data > libraries ----
  tar_target(libraries_file, "Libraries.tsv", format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file)),
  tar_target(libraries_col_data, make_col_data(list(libraries, samples))),

  # assay data (counts) ----
  tar_target(counts_file, find_counts_file(data_dir_name), format = "file"),
  tar_target(counts, readr::read_tsv(counts_file) |> tidy_counts()),
  tar_target(assay_data, make_assay_data(counts)),

  # summary data from previous steps ----
  tar_target(libraries_summary_file, find_libraries_summary_file(data_dir_name), format = "file"),
  tar_target(libraries_summary, readr::read_tsv(libraries_summary_file) |> tidy_libraries_summary()),
  tar_target(libraries_summary_rows, make_previous_summary_rows(libraries_summary, libraries_summary_file, base_provenance)),

  # row data (taxonomy) ----
  tar_target(taxonomy_file, find_taxonomy_file(data_dir_name, config |> pluck("taxonomy")), format = "file"),
  tar_target(taxonomy, readr::read_tsv(taxonomy_file) |> tidy_taxonomy()),
  tar_target(ranks, detect_taxonomy_ranks(taxonomy)),
  tar_target(features_info_file, find_features_info_file(data_dir_name), format = "file"),
  tar_target(features_info, readr::read_tsv(features_info_file) |> tidy_features_info()),
  tar_target(row_data, make_row_data(taxonomy, features_info)),

  # SummarizedExperiment > libraries > raw ----
  tar_target(se_libs_raw_provenance, modifyList(base_provenance, list(resolution = "libraries", state = "raw"))),
  tar_target(se_libs_raw, make_se(assay_data, libraries_col_data, row_data, ranks, se_libs_raw_provenance) |> add_decontam(negative_controls)),

  # SummarizedExperiment > samples > raw ----
  tar_target(se_raw, merge_cols(se_libs_raw, samples |> names() |> head(1L), samples |> names(), list(resolution = "samples"))),

  # SummarizedExperiment > samples > refined ----
  # filter features:
  tar_target(
    se_refined,
    se_raw |>
      filter_contaminants(decontam_threshold) |>
      keep_desirable_features(desirables) |>
      filter_undesirable_features(undesirables) |>
      update_provenance(se_raw, list(state = "refined")) |>
      tidy()
  ),
  tar_target(se_refined_filtered_features_table, filtered_features_table(se_refined)),
  tar_target(
    se_refined_filtered_features_file,
    write_tsv(
      se_refined_filtered_features_table,
      fs::path(results_dir_name, stringr::str_c(file_prefix, "filtered_features", sep = "_"), ext = "tsv"),
      na = "explicit"
    ),
    format = "file"
  ),
  # filter samples:
  tar_target(
    se_deep,
    se_refined |>
      filter_samples_by_sum(yield_min, yield_max) |>
      tidy()
  ),

  # SummarizedExperiment > all ----
  tar_target(se, list(se_libs_raw, se_raw, se_refined, se_deep)),
  tar_target(se_flat_file, export_flattened(se, results_dir_name), format = "file", pattern = map(se)),
  tar_target(ps, as_phyloseq(se), pattern = map(se)),

  ## SummarizedExperiment > all > summary ----
  tar_target(se_summary_rows, summary_as_row(se), pattern = map(se)),
  tar_target(
    se_summary,
    se_summary_rows |>
      dplyr::bind_rows() |>
      dplyr::bind_rows(libraries_summary_rows, se_part = _) |>
      dplyr::relocate(
        sample, feature, samples, features, total_counts, min_sample_counts, max_sample_counts, median_sample_counts,
        .after = last_col()
      )
  ),
  tar_target(
    se_summary_file,
    write_tsv(se_summary, fs::path(results_dir_name, stringr::str_c(file_prefix, "summary", sep = "_"), ext = "tsv")),
    format = "file"
  ),

  # beta diversity ----
  tar_target(ps_distance, calulcate_distance(ps), pattern = map(ps)),
  ## tests ----
  tar_target(permanova,
    test_distance(ps_distance, variable_of_interest),
    pattern = cross(ps_distance, variable_of_interest)
  ),
  tar_target(permanovas, permanova |> bind_rows()),
  tar_target(permanovas_file,
    write_tsv(permanovas, fs::path(results_dir_name, stringr::str_c(file_prefix, "tests", sep = "_"), ext = "tsv")),
    format = "file"
  ),
  ## ordination ----
  tar_target(ps_ordination, calulcate_ordination(ps_distance), pattern = map(ps_distance)),
  tar_target(ordination_plot_raw,
    plot_ordination(ps_ordination, variable_of_interest, sample_label_from, limits, theme),
    pattern = cross(ps_ordination, variable_of_interest)
  ),
  tar_target(ordination_plot,
    plot_ordination_with_tests(ordination_plot_raw, permanova),
    pattern = map(ordination_plot_raw, permanova)
  ),
  tar_target(ordination_plot_file, save_plot(ordination_plot, plots_dir_name), format = "file", pattern = map(ordination_plot)),

  # summary report ----
  tar_target(
    summary_report,
    make_summary_report(
      base_provenance,
      pipeline_version,
      list(
        samples = samples_file,
        libraries = libraries_file,
        counts = counts_file,
        `libraries summary` = libraries_summary_file,
        taxonomy = taxonomy_file
      ),
      list(
        desirables = desirables,
        undesirables = undesirables,
        yield_min = yield_min
      ),
      config,
      se_summary
    )
  ),
  tar_target(
    summary_report_file,
    write_text(summary_report, fs::path(results_dir_name, stringr::str_c(file_prefix, "summary", sep = "_"), ext = "md")),
    format = "file"
  )
)
