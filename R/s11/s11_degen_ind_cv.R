# s11_degen_ind_cv.R -- S11 (D-90, D-89 supplement): the per-replicate usable indicators of S10 "Degeneri"
# (results_v4/S10/degen/ind, one row per usable replicate) with the D-90 no-result datasets removed from the MMM with
# the convergence check (T5_mmm_cc: all four mechanisms; T34_mmm: main grid). The other files are symlinked unchanged.
# Read-only on results_v4/S10. Output: results_v4/S11/degen/ind/, tables/V8_S11_degen_ind_cv.csv.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L); setwd("/home/user/MMM")
source("R/s11/s11_cv.R")
SRC <- "results_v4/S10/degen/ind"; OUT <- "results_v4/S11/degen/ind"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
R <- list()
for (f in list.files(SRC)) {
  if (f %in% c("T5_mmm_cc.rds", "T34_mmm.rds")) {
    Y <- readRDS(file.path(SRC, f)); n0 <- nrow(Y); Y[, keep := TRUE]
    s11_drop_cv(Y, "keep", rows = rep(TRUE, nrow(Y)), what = f)   # D-90
    Y <- Y[keep == TRUE][, keep := NULL]; saveRDS(Y, file.path(OUT, f))
    R[[f]] <- data.table(file = f, rows_before = n0, rows_after = nrow(Y), removed = n0 - nrow(Y))
  } else if (!file.exists(file.path(OUT, f))) file.symlink(normalizePath(file.path(SRC, f)), file.path(OUT, f))
}
X <- rbindlist(R); print(X); fwrite(X, "tables/V8_S11_degen_ind_cv.csv")
