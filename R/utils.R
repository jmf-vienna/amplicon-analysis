fortify <- function(x) {
  withr::local_collate("en_US.UTF-8")

  x |>
    as.factor() |>
    forcats::fct_na_value_to_level("N/A") |>
    droplevels()
}

first_name <- function(x, ...) {
  x |>
    names() |>
    dplyr::first(...)
}

first_id_name <- function(x) {
  x |>
    names() |>
    stringr::str_subset("[^A-Z]ID$") |>
    dplyr::first()
}

as_full_tibble <- function(x, var = "rowname") {
  x |>
    as.data.frame() |>
    tibble::rownames_to_column(var) |>
    tibble::as_tibble()
}

smart_bind_rows <- function(x) {
  widest_names <-
    x |>
    purrr::map(names) |>
    purrr::map(length) |>
    unlist() |>
    sort() |>
    tail(1L)

  ordered_names <-
    x |>
    purrr::chuck(names(widest_names)) |>
    names()

  x |>
    dplyr::bind_rows() |>
    dplyr::relocate(tidyselect::all_of(ordered_names))
}

smart_arrange <- function(x, cols = c("tool", "resolution", "state", "phase")) {
  x |>
    dplyr::mutate(dplyr::across(any_of(cols), forcats::fct_inorder)) |>
    dplyr::arrange(dplyr::across(any_of(cols)))
}
