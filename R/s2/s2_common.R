# s2_common.R -- shared setup for the S2 scripts.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep),
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressMessages({
  library(data.table); library(MMLB); library(lme4); library(geepack); library(geesmv)
  library(dplyr); library(tibble); library(MASS); library(fastGHQuad)
})
REPO   <- path.expand("~/MMM")
RUNDIR <- file.path(REPO, ".rundir")
OUTDIR <- file.path(REPO, "results_v4", "S2"); dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TABDIR <- file.path(REPO, "tables")
source(file.path(REPO, "R/protocol/dgm_common.R"))
N_LEVELS <- c(10, 15, 20, seq(25, 145, 10)); T_LEVELS <- c(2, 3, 4, 6, 8, 10)
RNG_PROD <- c("L'Ecuyer-CMRG", "Inversion", "Rejection")   # D-07
set_rep_seed <- function(model, i) { RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(seed_of(model, i)) }
mcse_p <- function(p, n) sqrt(p * (1 - p) / n)
