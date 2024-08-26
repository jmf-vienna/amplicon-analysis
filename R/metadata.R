variable_info <- function(x) {
  counts <-
    vctrs::vec_count(x, sort = "none") |>
    dplyr::filter(jmf::a(key))

  all <- counts |>
    dplyr::pull(count)
  gte2 <- counts |>
    dplyr::filter(count >= 2L) |>
    dplyr::pull(count)
  gte3 <- counts |>
    dplyr::filter(count >= 3L) |>
    dplyr::pull(count)

  list(
    unique = all |> length(),
    gte2 = list(
      unique = gte2 |> length()
    ),
    gte3 = list(
      unique = gte3 |> length()
    )
  )
}
