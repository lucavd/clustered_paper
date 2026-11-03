# s6_mmm_iterlim.R -- intervention test of the mechanism (post hoc). Hypothesis: with a latent-SD slope
# near the boundary, nlm spends its 100-iteration budget moving the slope towards -Inf and leaves
# beta short of the maximum; the larger N, the sharper the likelihood in beta, the larger the error
# relative to the SE. Test: refit the NON-converged study fits of M3-M5, T >= 8, N = 55 and 145
# (from results_v4/S6/mmm_boundary_all.rds) with iter.lim = 2000, everything else as in the study.
# Same existing datasets (regenerated with their seeds). Output: results_v4/S6/mmm_iterlim_<dgp>.rds
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
FORM <- list(`3` = list(m = binY ~ time * trt, lv = ~time, t = ~1), `4` = list(m = binY ~ time * trt, lv = ~time, t = ~trt),
             `5` = list(m = binY ~ time * trt + dose, lv = ~time + dose, t = ~trt + steroids))
X <- readRDS(file.path(REPO, "results_v4/S6/mmm_boundary_all.rds"))
J <- X[dgp == DGP & useA == TRUE & nonconv == TRUE & model %in% 3:5 & ntime >= 8 & N %in% c(55L, 145L), .(model, ntime, N, i)]
J <- J[order(-N * ntime)]
one <- function(r) {
  cc <- J[r]
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 1000L * cc$model + cc$i)
  d <- as.data.frame(make_all_data(cc$N, cc$ntime, cc$i, cc$model, "manual")); F <- FORM[[as.character(cc$model)]]
  t0 <- proc.time()[[3]]
  x <- tryCatch(suppressWarnings(mm(F$m, lv.formula = F$lv, t.formula = F$t, data = d, id = id, q = 20, iter.lim = 2000)), error = function(e) NULL)
  if (is.null(x)) return(data.table(cc, L_ok = FALSE))
  k <- which(names(x$beta) == "time:trt")
  data.table(cc, L_ok = TRUE, L_est = unname(x$beta[k]), L_se_mod = sqrt(x$mod.cov[k, k]), L_code = unname(x$control["convergence_code"]),
             L_niter = unname(x$control["n_iter"]), L_alpha = paste(sprintf("%.6g", x$alpha), collapse = ";"), L_loglik = as.numeric(x$logLik), L_sec = proc.time()[[3]] - t0)
}
OUT <- rbindlist(parallel::mclapply(seq_len(nrow(J)), one, mc.cores = CORES, mc.preschedule = FALSE), fill = TRUE)
OUT[, dgp := DGP]
saveRDS(OUT, file.path(REPO, sprintf("results_v4/S6/mmm_iterlim_%s.rds", DGP)))
cat(DGP, nrow(OUT), "refits\n")
