# =============================================================================
# s4_common.R -- shared setup for S4 (re-aggregation with MCSE).
# Usability classes (U0/U1/U2) are coded exactly as in R/s2/s2_ab_usable_missing.R.
# =============================================================================
REPO <- path.expand("~/MMM")
LIB  <- file.path(REPO, "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep),
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressMessages(library(data.table))
TABDIR <- file.path(REPO, "tables")
OUT_S4 <- file.path(REPO, "results_v4/S4"); dir.create(OUT_S4, recursive = TRUE, showWarnings = FALSE)
LOG_S4 <- file.path(REPO, "reports/morris/S4_logs"); dir.create(LOG_S4, recursive = TRUE, showWarnings = FALSE)
# Variant runs (D-32): MMM_VARIANT=<tag> redirects every input/output path of S4-S6 to a parallel tree,
# so a reanalysed row set goes through the same scripts without touching the study outputs. Unset = unchanged.
VARIANT <- Sys.getenv("MMM_VARIANT", "")
if (nzchar(VARIANT)) {
  TABDIR <- file.path(REPO, "tables", paste0("variant_", VARIANT))
  OUT_S4 <- file.path(REPO, "results_v4", paste0("variant_", VARIANT), "S4")
  LOG_S4 <- file.path(REPO, "reports/morris", paste0("variant_", VARIANT), "S4_logs")
  for (p in c(TABDIR, OUT_S4, LOG_S4)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
WALD  <- c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink")
LRT   <- "glmer_lrt"
PRIM  <- c(WALD, LRT)                       # the 8 primary methods (intersection set)
SENS  <- "glmer_bobyqa"                     # D-16 sensitivity, from glmer_lrt rows
METHODS <- c(PRIM, SENS)
MARGINAL_METHODS <- c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "shrink")
COND_METHODS     <- c("glmer", "glmer_bobyqa")
Z975 <- qnorm(0.975)
CELL <- c("dgp", "scen", "N", "ntime", "label")
DP   <- c("dgp", "scen", "N", "ntime")
mcse_p <- function(p, n) sqrt(p * (1 - p) / n)

# U0 filter of agg_all5.R (Wald rows need a usable SE, LRT rows only a p-value)
u0_flag <- function(method, converged, estimate, se, p) {
  w <- converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p)
  l <- converged & is.finite(p)
  out <- ifelse(method == LRT, l, w); out[is.na(out)] <- FALSE; out
}
# warning classes of S2 (R/s2/s2_ab_usable_missing.R, CLS and 'primary')
CLS <- c(w_grad = "max\\|grad\\|", w_hess = "Hessian|unable to evaluate scaled gradient",
         w_unident = "nearly unidentifiable", w_optim = "failure to converge in|convergence code|bobyqa|Nelder_Mead",
         w_glm01 = "fitted probabilities numerically 0 or 1", w_glmnc = "algorithm did not converge",
         w_nan = "NaNs produced", w_nared = "NA/NaN replaced by maximum positive value")
conv_warn_flag <- function(status) {
  s <- ifelse(is.na(status), "", status)
  lme4c <- grepl(CLS[["w_grad"]], s) | grepl(CLS[["w_hess"]], s) | grepl(CLS[["w_unident"]], s) | grepl(CLS[["w_optim"]], s)
  glmnc <- grepl(CLS[["w_glmnc"]], s)
  !is.na(status) & (lme4c | glmnc)       # primary class lme4_convergence or glm_not_converged
}
add_usability <- function(A) {
  A[, U0 := u0_flag(method, converged, estimate, se, p)]
  A[, U1 := U0 & is.na(status)]
  A[, U2 := U0 & !conv_warn_flag(status)]
  invisible(A)
}
