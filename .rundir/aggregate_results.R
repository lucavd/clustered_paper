# =============================================================================
# aggregate_results.R
# Collapse the per-chunk replicate files into the paper's result tables.
#
# Reads ~/MMM/results_full/*.rds (one file per cell x 100-replicate chunk) and
# produces, for every (dgp, scenario, N, T, method) cell:
#   - ADEMP performance measures with Monte Carlo SEs   (performance_metrics.R)
#   - convergence / exclusion accounting
#   - both scoring rules: each method against its OWN estimand, and all methods
#     against the COMMON marginal truth
#
# Designed to run on a PARTIAL result set: cells with too few replicates are
# reported with their actual n_rep rather than silently averaged in, so the
# tables can be inspected while the simulation is still running.
#
# Usage:  Rscript aggregate_results.R [results_dir] [out_dir]
# =============================================================================

suppressMessages({library(dplyr); library(tidyr)})

`%||%` <- function(a, b) if (is.null(a)) b else a

args      <- commandArgs(trailingOnly = TRUE)
RESULTDIR <- if (length(args) >= 1) args[1] else path.expand("~/MMM/results_full")
OUTDIR    <- if (length(args) >= 2) args[2] else path.expand("~/MMM/tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# performance_metrics.R supplies summarise_ademp() / summarise_convergence().
# Look next to this script and in the working directory; do NOT inspect the call
# stack -- sys.frame() errors out when the script is run directly by Rscript
# rather than source()d.
if (!exists("summarise_ademp")) {
  here <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=",
            commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
  for (p in unique(c("performance_metrics.R",
                     file.path(here %||% ".", "performance_metrics.R"),
                     file.path("R", "performance_metrics.R")))) {
    if (!is.na(p) && file.exists(p)) { source(p); break }
  }
  if (!exists("summarise_ademp"))
    stop("performance_metrics.R not found next to aggregate_results.R")
}

# ---- load ------------------------------------------------------------------
files <- list.files(RESULTDIR, pattern = "[.]rds$", full.names = TRUE)
if (!length(files)) stop("no result chunks found in ", RESULTDIR)
cat(sprintf("reading %d chunk files from %s\n", length(files), RESULTDIR))

read_one <- function(f) {
  d <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d
}
raw <- bind_rows(lapply(files, read_one))
cat(sprintf("replicate rows: %s | cells present: %d\n",
            format(nrow(raw), big.mark = ","),
            nrow(distinct(raw, dgp, model, N, ntime, label, method))))

# ---- completeness map ------------------------------------------------------
# Which grid cells are finished, and which are still partial? Reported so no
# table is read as complete when it is not.
coverage <- raw %>%
  filter(term == "time:trt") %>%
  count(dgp, model, N, ntime, label, method, name = "n_rep") %>%
  mutate(complete = n_rep >= 1000)
write.csv(coverage, file.path(OUTDIR, "coverage_map.csv"), row.names = FALSE)
cat(sprintf("cells at full 1000 replicates: %d / %d (%.1f%%)\n",
            sum(coverage$complete), nrow(coverage),
            100 * mean(coverage$complete)))

# ---- ADEMP tables, both scoring rules --------------------------------------
# H1 cells (label == "mmlb") carry the effect; H0 cells give type I error.
for (sc in c("own", "marginal")) {
  tab <- summarise_ademp(raw, term_of_interest = "time:trt", scoring = sc)
  if (is.null(tab)) { cat("no rows for scoring =", sc, "\n"); next }
  f <- file.path(OUTDIR, sprintf("ademp_%s.csv", sc))
  write.csv(tab, f, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", basename(f), nrow(tab)))
}

# ---- type I error from the H0 arm ------------------------------------------
# Under H0 the tested coefficient is 0 while intercept and main effects match H1
# (see generation_model1.R), so the rejection rate is a like-for-like type I
# error and is directly comparable to the power computed from the H1 arm.
type1 <- raw %>%
  filter(term == "time:trt", label == "manual", converged, is.finite(p)) %>%
  group_by(dgp, model, N, ntime, method) %>%
  summarise(n_rep = n(),
            type1 = mean(p < 0.05),
            type1_mcse = sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / n()),
            .groups = "drop") %>%
  mutate(across(c(type1, type1_mcse), ~ round(100 * .x, 2)))
write.csv(type1, file.path(OUTDIR, "type1_error.csv"), row.names = FALSE)
cat(sprintf("wrote type1_error.csv (%d rows)\n", nrow(type1)))

# ---- convergence / exclusions ----------------------------------------------
conv <- summarise_convergence(raw, term_of_interest = "time:trt")
write.csv(conv, file.path(OUTDIR, "convergence.csv"), row.names = FALSE)

# geesmv timeouts are a finding in their own right: the corrections are not
# merely slow at small N, they fail to return at all on sparse data.
timeouts <- raw %>%
  filter(grepl("^gee_(md|kc|fg)$", method)) %>%
  group_by(dgp, N, ntime, method) %>%
  summarise(n = n(),
            pct_timeout = round(100 * mean(grepl("timeout", status %||% ""), na.rm = TRUE), 2),
            .groups = "drop")
write.csv(timeouts, file.path(OUTDIR, "geesmv_timeouts.csv"), row.names = FALSE)
cat(sprintf("wrote convergence.csv and geesmv_timeouts.csv\n"))

saveRDS(raw, file.path(OUTDIR, "raw_all.rds"), compress = "xz")
cat("done\n")

#' Pool per-cell ADEMP measures across a grid, dropping cells that carry no
#' estimate at all.
#'
#' Cells where a method could not be fitted (e.g. a random-slope GLMM at two
#' measurement occasions -- structurally unidentifiable) return NA measures. A
#' plain mean() over cells then propagates NA and wipes out the whole row, which
#' is how the first GenBinaryY summary reported NA for the GLMM despite 419 of
#' 423 cells being perfectly fine. Drop them explicitly AND report how many
#' cells contributed, so a shrinking denominator is visible rather than silent.
pool_ademp <- function(A, by = c("method", "estimand")) {
  A %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::summarise(
      cells_total    = dplyr::n(),
      cells_with_est = sum(!is.na(bias)),
      truth      = dplyr::first(truth),
      bias       = mean(bias, na.rm = TRUE),
      rmse       = mean(rmse, na.rm = TRUE),
      emp_se     = mean(emp_se, na.rm = TRUE),
      coverage   = 100 * mean(coverage, na.rm = TRUE),
      se_ratio   = mean(se_ratio, na.rm = TRUE),
      mape       = mean(mape, na.rm = TRUE),
      pct_separated = mean(pct_separated, na.rm = TRUE),
      .groups = "drop"
    )
}
