make_se_subsets <- function(se, subset) {
  if (is.null(subset)) {
    return(se)
  }

  loadNamespace(class(se))

  variable <- subset |> purrr::chuck("variable")
  values <- subset |> purrr::chuck("values")

  name <- subset |> purrr::pluck("name", .default = str_c(variable, ": ", values |> str_c(collapse = ", ")))

  # in case of multiple subsetting
  name <-
    se |>
    get_provenance() |>
    purrr::pluck("subset") |>
    c(name) |>
    str_c(collapse = " & ")

  se[, se[[variable]] %in% values] |>
    update_provenance(new = list(subset = name)) |>
    tidy()
}
