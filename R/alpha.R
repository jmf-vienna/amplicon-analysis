add_alpha_diversity <- function(se, alpha_diversity_indexes, threshold, rarefaction_rounds = 10L) {
  min <- min_col_sum(se)
  if (min >= threshold) {
    se <-
      se |>
      mia::addAlpha(
        index = alpha_diversity_indexes,
        name = str_c(".alpha_diversity_", alpha_diversity_indexes)
      ) |>
      mia::addAlpha(
        index = alpha_diversity_indexes,
        name = str_c(".alpha_diversity_at_", min, "_", alpha_diversity_indexes),
        sample = min,
        niter = rarefaction_rounds
      )
  }
  se
}
