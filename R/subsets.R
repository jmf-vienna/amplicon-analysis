make_se_subsets <- function(se, subset) {
  if (is.null(subset)) {
    return(se)
  }

  loadNamespace(class(se))

  variable <- subset |> chuck("variable")
  values <- subset |> chuck("values")
  name <- subset |> pluck("name", .default = str_c(variable, ": ", values |> str_c(collapse = ", ")))

  se[, se[[variable]] %in% values] |>
    update_provenance(new = list(subset = name)) |>
    tidy()
}
