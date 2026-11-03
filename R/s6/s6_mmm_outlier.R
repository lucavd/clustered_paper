# s6_mmm_outlier.R -- what distinguishes the extreme MMM estimates in M3 (GenBinaryY, H0, N = 145, T = 10)?
# Same 200 existing replicates as s6_mmm_se_probe.R (regenerated with their own seeds, no new data).
# Refits mm() and keeps the convergence code, the iterations, the latent-variable (alpha) and
# transition (gamma) estimates; compares replicates with |est - median| > 4 MAD against the rest,
# and refits every replicate also with a larger quadrature (q = 50).
# Output: tables/V4_S6_mmm_outlier.csv
source(path.expand("~/MMM/R/s6/s6_common.R"))
suppressMessages({ library(MMLB); library(dplyr) })
RD <- file.path(REPO, "R/v4/rundir")
for (f in c("generation_model3.R", "generate_data.R", "make_all_data.R")) source(file.path(RD, f))
P <- fread(file.path(TABDIR, "V4_S6_mmm_se_probe.csv"))[model == 3 & is.finite(se_mod)]
m <- median(P$est); s <- mad(P$est); P[, extreme := abs(est - m) > 4 * s]
one <- function(i) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 3000L + i)
  d <- as.data.frame(make_all_data(145L, 10L, i, 3L, "manual"))
  f <- function(q) tryCatch(mm(binY ~ time * trt, lv.formula = ~time, t.formula = ~1, data = d, id = id, q = q), error = function(e) NULL)
  a <- f(20); b <- f(50)
  g <- function(x, nm) if (is.null(x)) NULL else {
    k <- which(names(x$beta) == "time:trt")
    data.table(i = i, fit = nm, est = unname(x$beta[k]), se = sqrt(x$mod.cov[k, k]),
               code = unname(x$control["convergence_code"]), niter = unname(x$control["n_iter"]),
               gamma0 = unname(x$alpha[1]), lsig0 = unname(x$alpha[2]), lsig_time = unname(x$alpha[3]), loglik = as.numeric(x$logLik))
  }
  rbind(g(a, "q20"), g(b, "q50"))
}
OUTF <- file.path(TABDIR, "V4_S6_mmm_outlier.csv")
if (Sys.getenv("S6_SUMMARY_ONLY") == "1" && file.exists(OUTF)) R <- fread(OUTF) else {
  R <- rbindlist(parallel::mclapply(P$i, one, mc.cores = 30, mc.preschedule = FALSE), fill = TRUE)
  R <- merge(R, P[, .(i, extreme)], by = "i")
  fwrite(R, OUTF)
}
print(R[, .(n = .N, code4 = sum(code == 4), code_other = sum(code != 0 & code != 4), med_est = median(est), sd_est = sd(est), med_lsig_time = median(lsig_time),
           n_lsig_time_lt_m3 = sum(lsig_time < -3), med_niter = as.numeric(median(niter))), by = .(fit, extreme)], digits = 3)
print(R[fit == "q20", .(n = .N, sd_est = sd(est), med_se = median(se), rej = 100 * mean(abs(est / se) > qnorm(.975), na.rm = TRUE)), by = .(code)], digits = 3)
W <- dcast(R, i + extreme ~ fit, value.var = "est"); print(W[, .(sd_q20 = sd(q20), sd_q50 = sd(q50, na.rm = TRUE), maxdiff = max(abs(q20 - q50), na.rm = TRUE)), by = extreme], digits = 3)
print(R[extreme == TRUE & fit == "q20", .(i, est, se, code, niter, gamma0, lsig0, lsig_time)])
