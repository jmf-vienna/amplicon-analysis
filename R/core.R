get_pipeline_version <- function() {
  last_commit <- gert::git_log(max = 1L, repo = targets::tar_config_get("script"))
  hash <- last_commit[["commit"]] |> stringr::str_sub(1, 7)

  glue::glue("0.1.0-{hash}")
}
