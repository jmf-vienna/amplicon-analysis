as_tibble <- function(x, var = "rowname") {
  x |>
    as.data.frame() |>
    tibble::rownames_to_column(var) |>
    tibble::as_tibble()
}

bind_rows <- function(x) {
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
