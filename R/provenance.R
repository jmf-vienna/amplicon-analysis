get_provenance <- function(x) {
  xfun::attr(x, "provenance") |> as.list()
}

set_provenance <- function(x, provenance) {
  attr(x, "provenance") <- provenance
  invisible(x)
}

update_provenance <- function(x, source, new = list()) {
  get_provenance(x) |>
    modifyList(get_provenance(source)) |>
    modifyList(new) |>
    set_provenance(x, provenance = _)
}

as_file_name <- function(x) {
  x |>
    get_provenance() |>
    purrr::map(\(x) {
      stringr::str_c(names(x), x, sep = "_", collapse = "_")
    }) |>
    stringr::str_c(collapse = "_") |>
    stringr::str_replace_all(" ", "_")
}

as_title <- function(x) {
  stringr::str_c(names(x), x, sep = ": ", collapse = " | ")
}

plot_titles <- function(plot, data, n = 4) {
  plot + ggplot2::labs(
    title = data |> get_provenance() |> head(n) |> as_title() |> stringr::str_remove("^[a-z]+: "),
    subtitle = data |> get_provenance() |> tail(-n) |> as_title()
  )
}
