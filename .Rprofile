source("renv/activate.R")

options(tidyverse.quiet = TRUE)

if (
  interactive() &&
  requireNamespace("usethis", quietly = TRUE)
) {
  .attach <- function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      pkg |>
        library(character.only = TRUE) |>
        suppressWarnings() |>
        suppressPackageStartupMessages()
      usethis::ui_done("Package {ui_value(pkg)} attached.")
    }
  }
  usethis::ui_todo(
    "Attaching development supporting packages...
    (on interactive section only!)"
  )
  .attach("usethis")
  .attach("targets")
  .attach("tarchetypes")

  .run <- function(...) {
    source(here::here("dev/run.R")) # check and make pipeline.
    .run(...)
  }
  ui_info(
    "Execute {ui_code('.run()')} to make the pipeline interactively."
  )

  .background_run <- function(...) {
    stopifnot(requireNamespace("rstudioapi"))
    rstudioapi::jobRunScript(
      here::here("dev/background_run.R"),
      workingDir = here::here()
    )
  }
  ui_info(paste0(
    "Execute {ui_code('.background_run()')} to make the pipeline ",
    "as a background job in RStudio."
  ))
}
