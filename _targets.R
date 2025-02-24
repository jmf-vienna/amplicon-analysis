library(targets)

jmf::quiet()
options(warn = 2L)
targets::tar_option_set(
  packages = c("cli", "fs", "purrr", "readr", "rlang", "stringr", "vctrs"),
  format = "qs",
  iteration = "list",
  controller = crew::crew_controller_local(workers = 2L)
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

if (fs::file_exists("customize.R")) {
  source("customize.R")
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
  tar_target(rd_dir_name, config |> pluck("path", "rd", .default = "rd")),
  tar_target(file_prefix, base_provenance |> as_file_name()),

  ## config > refinement ----
  tar_target(negative_controls, config |> pluck("filter", "negative controls", .default = character())),
  tar_target(desirables, config |> pluck("filter", "desirable", .default = list())),
  tar_target(undesirables, config |> pluck("filter", "undesirable", .default = list())),
  tar_target(yield_min, config |> pluck("filter", "yield", "min", .default = 0L)),
  tar_target(yield_max, config |> pluck("filter", "yield", "max", .default = Inf)),
  tar_target(pass_libraries_yield_min, config |> pluck("filter", "pass", "yield", "min", .default = 0L)),
  tar_target(decontam_threshold, config |> pluck("filter", "decontam", .default = 0.0)),
  tar_target(length_min, config |> pluck("filter", "length", "min", .default = 0L)),
  tar_target(length_max, config |> pluck("filter", "length", "max", .default = Inf)),
  tar_target(failed_samples, config |> pluck("filter", "failed samples", .default = character())),

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
  tar_target(samples, readr::read_tsv(samples_file) |> tidy_samples()),
  tar_target(samples_col_data, make_col_data(list(samples))),

  # column data > sublibraries ----
  tar_target(sublibraries_file, "Sublibraries.tsv", format = "file"),
  tar_target(sublibraries, readr::read_tsv(sublibraries_file) |> tidy_sublibraries()),

  # column data > libraries ----
  tar_target(libraries_file, "Libraries.tsv", format = "file"),
  tar_target(libraries, readr::read_tsv(libraries_file) |> tidy_libraries()),
  tar_target(libraries_col_data, make_col_data(list(sublibraries, libraries, samples))),

  # assay data (counts) ----
  tar_target(counts_file, find_counts_file(data_dir_name), format = "file"),
  tar_target(counts, readr::read_tsv(counts_file) |> tidy_counts()),
  tar_target(assay_data, make_assay_data(counts)),

  # metrics from previous steps ----
  tar_target(prior_library_metrics_file, find_prior_metrics_file(data_dir_name), format = "file"),
  tar_target(
    prior_library_metrics,
    prior_library_metrics_file |>
      read_tsv() |>
      tidy_prior_metrics() |>
      dplyr::bind_cols(base_provenance |> tibble::as_tibble(), second_argument = _)
  ),

  # row data ----
  ## features info ----
  tar_target(features_info_file, find_features_info_file(data_dir_name), format = "file"),
  tar_target(features_info, readr::read_tsv(features_info_file) |> tidy_features_info()),
  ## taxonomy ----
  tar_target(taxonomy_file, find_taxonomy_file(data_dir_name, config |> pluck("taxonomy")), format = "file"),
  tar_target(taxonomy, readr::read_tsv(taxonomy_file, guess_max = Inf) |> tidy_taxonomy() |> taxonomy_fallback(features_info)),
  tar_target(ranks, detect_taxonomy_ranks(taxonomy)),
  ## merge ----
  tar_target(row_data, make_row_data(taxonomy, features_info)),

  # variable names ----
  tar_target(library_id_var, libraries |> first_name()),
  tar_target(biosample_id_var, samples |> first_name()),
  tar_target(feature_id_var, ranks |> dplyr::last()),

  # SEs ----
  ## library SEs ----
  ### raw ----
  tar_target(se_libs_raw_provenance, modifyList(base_provenance, list(resolution = "libraries", state = "raw"))),
  tar_target(se_libs_rawer, make_se(assay_data, libraries_col_data, row_data, ranks, se_libs_raw_provenance)),
  tar_target(failed_libraries, get_failed_libraries(se_libs_rawer, negative_controls, pass_libraries_yield_min, failed_samples)),
  tar_target(se_libs_raw, se_libs_rawer |> add_decontam(negative_controls, failed_libraries)),

  ### clean (filter features) ----
  tar_target(
    se_libs_clean,
    se_libs_raw |>
      filter_by_length(length_min, length_max) |>
      filter_contaminants(decontam_threshold) |>
      keep_desirable_features(desirables) |>
      filter_undesirable_features(undesirables) |>
      update_provenance(new = list(state = "clean")) |>
      tidy()
  ),

  ### filtered features table ----
  tar_target(filtered_features_table, make_filtered_features_table(se_libs_clean)),
  tar_target(
    filtered_features_file,
    write_tsv(
      filtered_features_table,
      fs::path(results_dir_name, stringr::str_c(file_prefix, "filtered_features", sep = "_"), ext = "tsv"),
      na = "explicit"
    ),
    format = "file"
  ),

  ### pass (filter failed libraries) ----
  tar_target(
    se_libs_pass,
    se_libs_clean |>
      remove_cols(failed_libraries) |>
      update_provenance(new = list(state = "pass")) |>
      tidy()
  ),

  ## sample SEs ----
  ### raw (not used later on) ----
  tar_target(se_raw, merge_cols(
    se_libs_raw,
    biosample_id_var,
    samples |> names(),
    list(resolution = "samples")
  )),

  ### refined ----
  tar_target(se_refined, merge_cols(
    se_libs_pass,
    biosample_id_var,
    samples |> names(),
    list(resolution = "samples", state = "refined")
  )),

  ### deep (filter samples) ----
  tar_target(
    se_deep,
    se_refined |>
      filter_samples_by_sum(yield_min, yield_max) |>
      tidy()
  ),

  ## filtered samples table ----
  tar_target(filtered_samples_table, make_filtered_samples_table(list(
    list(se_libs_raw, se_libs_clean),
    list(se_libs_clean, se_libs_pass),
    list(se_refined, se_deep)
  ))),
  tar_target(
    filtered_samples_table_file,
    write_tsv(
      filtered_samples_table,
      fs::path(results_dir_name, stringr::str_c(file_prefix, "filtered_samples", sep = "_"), ext = "tsv")
    ),
    format = "file"
  ),

  ## all SEs ----
  tar_target(lib_se, list(se_libs_raw, se_libs_clean)),
  tar_target(sam_se, list(se_raw, se_refined, se_deep)),
  tar_target(se, list_c(list(lib_se, sam_se))),
  tar_target(se_file, export_se(se, rd_dir_name), format = "file", pattern = map(se)),
  tar_target(se_flat_file, export_flattened(se, results_dir_name), format = "file", pattern = map(se)),
  tar_target(ps, as_phyloseq(se), pattern = map(se)),
  tar_target(ps_file, export_ps(ps, rd_dir_name), format = "file", pattern = map(ps)),

  # library metrics ----
  tar_target(se_library_metrics, make_library_metrics(lib_se), pattern = map(lib_se)),
  tar_target(
    library_metrics,
    se_library_metrics |>
      dplyr::bind_rows() |>
      dplyr::bind_rows(prior_library_metrics, se_part = _) |>
      dplyr::inner_join(
        dplyr::select(libraries, library_id = all_of(library_id_var), biosample_id = all_of(biosample_id_var)),
        by = "library_id"
      ) |>
      dplyr::relocate(count, .after = last_col()) |>
      dplyr::arrange(library_id)
  ),
  tar_target(
    library_metrics_file,
    write_tsv(library_metrics, fs::path(results_dir_name, stringr::str_c(file_prefix, "library_metrics", sep = "_"), ext = "tsv")),
    format = "file"
  ),

  # biosample metrics ----
  tar_target(se_biosample_metrics, make_biosample_metrics(sam_se), pattern = map(sam_se)),
  tar_target(
    biosample_metrics,
    se_biosample_metrics |>
      dplyr::bind_rows() |>
      dplyr::bind_rows(
        library_metrics |>
          dplyr::group_by(across(!c(library_id, count))) |>
          dplyr::summarise(count = sum(count), .groups = "drop"),
        second_argument = _
      ) |>
      dplyr::relocate(sample_id_var, biosample_id, libraries, count, .after = last_col()) |>
      dplyr::arrange(biosample_id)
  ),
  tar_target(
    biosample_metrics_file,
    write_tsv(biosample_metrics, fs::path(results_dir_name, stringr::str_c(file_prefix, "biosample_metrics", sep = "_"), ext = "tsv")),
    format = "file"
  ),

  # summary rows ----
  tar_target(se_summary_rows, summary_as_row(se), pattern = map(se)),
  tar_target(se_summary, se_summary_rows |> dplyr::bind_rows()),
  tar_target(metrics_summary, make_metrics_summary(library_metrics, biosample_metrics, se_summary)),
  tar_target(
    metrics_summary_file,
    write_tsv(metrics_summary, fs::path(results_dir_name, stringr::str_c(file_prefix, "summary", sep = "_"), ext = "tsv")),
    format = "file"
  ),

  # beta diversity ----
  tar_target(ps_distance, calulcate_distance(ps), pattern = map(ps)),
  ## tests ----
  tar_target(permanova,
    test_distance(ps_distance, variable_of_interest),
    pattern = cross(ps_distance, variable_of_interest)
  ),
  tar_target(permanovas, permanova |> bind_rows() |> finalize_tests_table()),
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
        `library metrics` = prior_library_metrics_file,
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
