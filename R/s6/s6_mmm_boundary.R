# =============================================================================
# s6_mmm_boundary.R -- S6, PI question on MMM with long series (post hoc diagnostic, pre-registered
# in docs/steps/S6.md before running). No new simulated data: the first R replicates of existing
# cells are regenerated with their own seeds (set.seed(20260805 + 1000 * model + i), L'Ecuyer-CMRG,
# as R/v4/rundir/run_full.R) and refitted.
#
# Cause under test. MMLB::mm exponentiates EACH latent-SD coefficient (LogLScoreCalc:
# sigma <- exp(params[sig]); sigma_it = Xlv %*% sigma), so with lv.formula = ~time
# sigma_t = exp(a0) + exp(a1) * t: the slope is forced positive and slope = 0 is the boundary
# a1 -> -Inf. All three generators have a latent-SD slope of about 0.01 per occasion, close to it.
#
# Per dataset: fit A = the study fit (formulas of R/v4/rundir/analysis_model*.R, q = 20, defaults),
# keeping est, model-based SE (mod.cov), robust SE (rob.cov), nlm code, iterations, all alpha;
# fit B = the same model without the latent-SD slope terms (lv.formula = ~1), the boundary-reduced fit.
# Usage: Rscript s6_mmm_boundary.R <dgp> <cores>   dgp in GenBinaryY | GLMM-native | GLMM-AR1
# Output: results_v4/S6/mmm_boundary_<dgp>.rds
# =============================================================================
args <- commandArgs(TRUE); DGP <- args[1]; CORES <- as.integer(args[2])
source(path.expand("~/MMM/R/s6/s6_common.R"))
suppressMessages({ library(MMLB); library(dplyr) })
RD <- file.path(REPO, "R/v4/rundir"); sp <- function(f) source(file.path(RD, f))
if (DGP == "GenBinaryY") {
  for (f in c(sprintf("generation_model%d.R", 1:5), "generate_data.R", "make_all_data.R")) sp(f)
} else if (DGP == "GLMM-native") {
  sp("gen_glmm_data.R"); for (f in c(sprintf("manual/generation_model%d.R", 1:5), "manual/generate_data.R", "manual/make_all_data.R")) sp(f)
} else {
  sp("gen_ar1_data.R"); for (f in c("ar1/generation_models_ar1.R", "manual/generate_data.R", "ar1/make_all_data.R")) sp(f)
}
FORM <- list(`1` = list(m = binY ~ time * trt, lv = ~1, t = ~1),
             `3` = list(m = binY ~ time * trt, lv = ~time, t = ~1),
             `4` = list(m = binY ~ time * trt, lv = ~time, t = ~trt),
             `5` = list(m = binY ~ time * trt + dose, lv = ~time + dose, t = ~trt + steroids))
R_REP <- as.integer(Sys.getenv("S6_RREP", "200"))
cells <- CJ(model = c(1L, 3L, 4L, 5L), ntime = c(2L, 4L, 6L, 8L, 10L), N = c(10L, 25L, 55L, 145L), i = seq_len(R_REP))
cells <- cells[order(-N * ntime)]                       # largest fits first, for load balance
fitone <- function(d, F, lv) {
  t0 <- proc.time()[[3]]
  x <- tryCatch(suppressWarnings(mm(F$m, lv.formula = lv, t.formula = F$t, data = d, id = id, q = 20)), error = function(e) NULL)
  if (is.null(x)) return(list(ok = FALSE, sec = proc.time()[[3]] - t0))
  k <- which(names(x$beta) == "time:trt")
  list(ok = TRUE, est = unname(x$beta[k]), se_mod = sqrt(x$mod.cov[k, k]), se_rob = sqrt(x$rob.cov[k, k]),
       code = unname(x$control["convergence_code"]), niter = unname(x$control["n_iter"]),
       alpha = paste(sprintf("%.6g", x$alpha), collapse = ";"), alpha_names = paste(names(x$alpha), collapse = ";"),
       loglik = as.numeric(x$logLik), sec = proc.time()[[3]] - t0)
}
one <- function(r) {
  cc <- cells[r]
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 1000L * cc$model + cc$i)
  d <- tryCatch(as.data.frame(make_all_data(cc$N, cc$ntime, cc$i, cc$model, "manual")), error = function(e) NULL)
  if (is.null(d)) return(data.table(cc, gen_ok = FALSE))
  F <- FORM[[as.character(cc$model)]]
  A <- fitone(d, F, F$lv)
  B <- if (cc$model == 1L) list(ok = NA) else fitone(d, F, ~1)
  names(A) <- paste0("A_", names(A)); names(B) <- paste0("B_", names(B))
  data.table(cc, gen_ok = TRUE, as.data.table(A), as.data.table(B))
}
t0 <- Sys.time()
OUT <- rbindlist(parallel::mclapply(seq_len(nrow(cells)), one, mc.cores = CORES, mc.preschedule = FALSE), fill = TRUE)
OUT[, dgp := DGP]
dir.create(file.path(REPO, "results_v4/S6"), showWarnings = FALSE)
saveRDS(OUT, file.path(REPO, sprintf("results_v4/S6/mmm_boundary_%s.rds", DGP)))
cat(DGP, nrow(OUT), "rows;", format(difftime(Sys.time(), t0, units = "mins")), "\n")
