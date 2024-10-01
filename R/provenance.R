get_provenance <- function(x) {
  xfun::attr(x, "provenance") |> as.list()
}

set_provenance <- function(x, provenance) {
  attr(x, "provenance") <- provenance
  invisible(x)
}

update_provenance <- function(x, source = NULL, new = list()) {
  get_provenance(x) |>
    modifyList(get_provenance(source)) |>
    modifyList(new) |>
    set_provenance(x, provenance = _)
}

provenance_as_tibble <- function(x) {
  x |>
    get_provenance() |>
    purrr::list_flatten(name_spec = "{outer} {inner}") |>
    tibble::as_tibble()
}

provenance_as_file_name <- function(x) {
  x |>
    get_provenance() |>
    as_file_name()
}

as_file_name <- function(x) {
  x |>
    purrr::map(\(x) {
      stringr::str_c(
        names(x),
        x |> stringr::str_remove(" \\(.+\\)$"), # nolint: nonportable_path_linter.
        sep = "_", collapse = "_"
      )
    }) |>
    stringr::str_flatten("_") |>
    stringr::str_replace_all(stringr::fixed("≤"), "lte") |>
    stringr::str_replace_all(stringr::fixed("≥"), "gte") |>
    force_valid_file_name()
}

as_title <- function(x) {
  x <-
    x |>
    unlist()

  n <-
    x |>
    names() |>
    stringr::str_replace_all(stringr::fixed("."), ": ")

  stringr::str_c(n, dplyr::if_else(n == "", "", ": "), x) |>
    stringr::str_replace("(≤|≥): ", "\\1") |>
    stringr::str_flatten(" | ")
}

plot_titles <- function(plot, title_n = 2L, title = NULL, subtitle = NULL, subtitles = NULL) {
  provenance <-
    plot |>
    get_provenance() |>
    purrr::list_assign(aesthetics = rlang::zap())

  # trim "project:" from the title
  names(provenance)[[1L]] <- ""

  plot + ggplot2::labs(
    title = provenance |>
      head(title_n) |>
      c(title) |>
      as_title(),
    subtitle = c(
      provenance |>
        tail(-title_n) |>
        c(subtitle) |>
        as_title(),
      subtitles
    ) |>
      stringr::str_flatten("\n")
  )
}
