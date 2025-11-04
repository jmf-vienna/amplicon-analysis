get_se_version <- function() {
  1L
}

get_pipeline_version <- function() {
  repo <- targets::tar_config_get("script")

  n_commits <- gert::git_log(max = 1e9L, repo = repo) |> nrow()

  last_commit <- gert::git_log(max = 1L, repo = repo)
  hash <- last_commit[["commit"]] |> stringr::str_sub(1L, 7L)

  glue::glue("0.{get_se_version()}.0-beta.{n_commits}+{hash}")
}
