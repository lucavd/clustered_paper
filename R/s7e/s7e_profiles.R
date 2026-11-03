# s7e_profiles.R -- S7e calibration of the copula dependence, BEFORE any method is fitted to copula data.
# Within-subject dependence profile of each mechanism of the study (M1, H0, T = 10, one large dataset of
# 20,000 subjects per mechanism, seed 20261023) and of candidate Gaussian copulas: Pearson correlation of
# standardised residuals (y - mu) / sqrt(mu (1 - mu)) by lag, mu = cell mean by (time, trt); plus the
# exact Pearson correlation of the copula from the bivariate normal (closed form, no simulation).
source(path.expand("~/MMM/R/s3/s3_common.R"))
suppressMessages(library(data.table))
CD <- file.path(REPO, "R/v4/rundir/copula")
lag_prof <- function(d) {
  d <- as.data.table(d)[, .(id, time, trt, y = binY)]
  d[, mu := mean(y), by = .(time, trt)]; d[, r := (y - mu) / sqrt(mu * (1 - mu))]
  p <- merge(d[, .(id, t1 = time, r1 = r)], d[, .(id, t2 = time, r2 = r)], by = "id", allow.cartesian = TRUE)[t2 > t1]
  p[, .(cor = mean(r1 * r2), npairs = .N), by = .(lag = t2 - t1)][order(lag)]
}
N <- 20000L; Tm <- 10L; out <- list()
for (dg in DGPS) {
  env <- load_env(dg, RUNDIR_V4); RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(20261023L)
  t0 <- Sys.time(); d <- env$make_all_data(N, Tm, 1L, 1L, "manual")
  x <- lag_prof(d); x[, `:=`(mech = dg, prev = mean(d$binY), secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))]; out[[dg]] <- x
}
cenv <- new.env(parent = globalenv())
for (f in c("gen_copula_data.R", "generation_models_copula.R")) sys.source(file.path(CD, f), envir = cenv)
cand <- CJ(floor = c(0, 0.15, 0.3), rho = c(0.5, 0.6, 0.7, 0.8, 0.9))
for (k in seq_len(nrow(cand))) {
  cenv$COPULA_RHO <- cand$rho[k]; cenv$COPULA_FLOOR <- cand$floor[k]
  RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(20261023L)
  d <- cenv$generate_data_model1(N, Tm, 1L, "manual")
  x <- lag_prof(d); x[, `:=`(mech = sprintf("copula floor=%.2f rho=%.1f", cand$floor[k], cand$rho[k]), prev = mean(d$binY), secs = NA_real_)]
  out[[length(out) + 1L]] <- x
}
P <- rbindlist(out)
fwrite(P, file.path(REPO, "tables/V4_S7e_dependence_profiles.csv"))
# summary: lag 1, lag 3, lag 9, mean over all pairs of a full 10-visit series (weights 10 - lag) and of a 4-visit series
S <- P[, .(prev = prev[1], lag1 = cor[lag == 1], lag3 = cor[lag == 3], lag9 = cor[lag == 9],
           mean_T10 = sum(cor * (10 - lag)) / sum(10 - lag), mean_T4 = sum((cor * (4 - lag))[lag <= 3]) / sum((4 - lag)[lag <= 3])), by = mech]
fwrite(S, file.path(REPO, "tables/V4_S7e_dependence_summary.csv"))
print(S, digits = 3)
