keep_desirables <- function(se, config) {
  se |>
    SummarizedExperiment::subset(
      Domain %in% config[["Domain"]]
    )
}

filter_undesirables <- function(se, config) {
  se |>
    SummarizedExperiment::subset(!(
      Order %in% config[["Order"]] |
        Family %in% config[["Family"]] |
        Sequence_ID %in% config[["Sequence_ID"]]
    ))
}
