set_provenance <- function(dest, src, new = list()) {
  provenance <- c(
    xfun::attr(src, "provenance"),
    new
  )

  attr(dest, "provenance") <- provenance
  dest
}

get_provenance <- function(x) {
  attr(x, "provenance")
}

collapse_provenance <- function(x, sep = "_") {
  x |>
    get_provenance() |>
    stringr::str_c(collapse = sep)
}
