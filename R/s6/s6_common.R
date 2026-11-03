# =============================================================================
# s6_common.R -- shared setup for S6 (tables and figures, M1c classification).
# Reuses the S5 engine unchanged (R/s5/s5_engine.R) and the S5 bootstrap aggregates.
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R")); source(file.path(REPO, "R/s5/s5_engine.R"))
OUT_S6 <- file.path(REPO, "results_v4/S6"); dir.create(OUT_S6, recursive = TRUE, showWarnings = FALSE)
LOG_S6 <- file.path(REPO, "reports/morris/S6_logs"); dir.create(LOG_S6, recursive = TRUE, showWarnings = FALSE)
FIGDIR <- file.path(REPO, "reports/morris/S6_figures"); dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
REPDIR <- file.path(REPO, "reports/morris")
if (nzchar(VARIANT)) {                                    # D-32 variant tree (see s4_common.R)
  REPDIR <- file.path(REPO, "reports/morris", paste0("variant_", VARIANT))
  OUT_S6 <- file.path(REPO, "results_v4", paste0("variant_", VARIANT), "S6")
  LOG_S6 <- file.path(REPDIR, "S6_logs"); FIGDIR <- file.path(REPDIR, "S6_figures")
  for (p in c(OUT_S6, LOG_S6, FIGDIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
TCODE <- function(t) fifelse(t <= 3, "T23", "T4p")
# M1c uses the model and the seed of M1: relabelled as M1 it goes through the main-grid engine unchanged
relabel_m1c <- function(cl) { x <- cl[scen == "M1c"]; x[, scen := "M1"]; x }
