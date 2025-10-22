cli::cli_alert("sourcing .Rprofile")
library(conflicted)

if (interactive()) {
  require(targets)

  require(cli)
  require(dplyr)
  require(forcats)
  require(fs)
  require(ggplot2)
  require(purrr)
  require(readr)
  require(rlang)
  require(stringr)
  require(tibble)
  require(tidyr)
  require(vctrs)

  tar_load_globals()
} else {
  cli::cli_abort("sourcing .Rprofile non-interactively")
}
