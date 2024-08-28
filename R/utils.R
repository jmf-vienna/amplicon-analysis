as_tibble <- function(x, var = "rowname") {
  x |>
    as.data.frame() |>
    tibble::rownames_to_column(var) |>
    tibble::as_tibble()
}
