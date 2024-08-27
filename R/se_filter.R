keep_desirables <- function(se, config) {
  se |>
    SummarizedExperiment::subset(
      SummarizedExperiment::rowData(se)$Domain %in% config[["Domain"]]
    )
}

filter_undesirables <- function(se, config) {
  se |>
    SummarizedExperiment::subset(!(
      SummarizedExperiment::rowData(se)$Order %in% config[["Order"]] |
        SummarizedExperiment::rowData(se)$Family %in% config[["Family"]] |
        SummarizedExperiment::rowData(se)$Sequence_ID %in% config[["Sequence_ID"]]
    ))
}
