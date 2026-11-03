# =============================================================================
# s2_e_corr.R -- S2(e), Q5: correlation of the time:trt estimate between arms
# (H0 vs H1, same mechanism) and between mechanisms (same arm), at equal
# (scenario, N, T, sim_id). Usable rows only (U0 = filter of agg_all5.R).
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
A <- readRDS(file.path(OUTDIR, "rows_timetrt.rds"))
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink")
A <- A[term == "time:trt" & method %in% WALD & converged & is.finite(estimate) & is.finite(se) &
       abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p), .(dgp, model, N, ntime, label, method, sim_id, estimate)]
fz <- function(r, n) { if (n < 4 || !is.finite(r)) return(c(NA_real_, NA_real_)); z <- atanh(r); tanh(z + c(-1, 1) * qnorm(0.975) / sqrt(n - 3)) }
corr_tab <- function(W, a, b) W[, {
  ok <- is.finite(get(a)) & is.finite(get(b)); n <- sum(ok)
  r <- if (n >= 4) cor(get(a)[ok], get(b)[ok]) else NA_real_; ci <- fz(r, n)
  .(n_pairs = n, rho = r, lo = ci[1], hi = ci[2])
}, by = setdiff(names(W), c(a, b, "sim_id"))]

# ---- Q5.1 between arms ---------------------------------------------------------
W1 <- dcast(A, dgp + model + N + ntime + method + sim_id ~ label, value.var = "estimate")
C1 <- corr_tab(W1, "manual", "mmlb")
fwrite(C1, file.path(TABDIR, "V4_S2_corr_arms.csv"))
k <- C1[method %in% c("mmm", "gee", "glmer") & n_pairs >= 100]
cat(sprintf("Q5.1 combos (mmm, gee, glmer; >= 100 pairs): %d; rho > 0.5: %d (%.1f%%); pass >= 95%%: %s\n",
            nrow(k), sum(k$rho > 0.5), 100 * mean(k$rho > 0.5), mean(k$rho > 0.5) >= 0.95))
cat(sprintf("Q5.1 lower CI > 0: %d of %d combos with >= 100 pairs (all methods: %d of %d); all combos with >= 4 pairs: %d of %d\n",
            sum(k$lo > 0), nrow(k), sum(C1[n_pairs >= 100]$lo > 0), nrow(C1[n_pairs >= 100]),
            sum(C1[n_pairs >= 4]$lo > 0), nrow(C1[n_pairs >= 4])))
cat("Q5.1 rho between arms, median [min, max] by dgp x method (>= 100 pairs):\n")
print(dcast(C1[n_pairs >= 100, .(v = sprintf("%.3f [%.3f, %.3f]", median(rho), min(rho), max(rho))), by = .(dgp, method)], method ~ dgp, value.var = "v"))
cat("Q5.1 combos with rho <= 0.5 (mmm, gee, glmer):\n"); print(k[rho <= 0.5][order(rho)][1:min(.N, 30)])

# ---- Q5.2 between mechanisms ---------------------------------------------------
A[, dgp := c("GenBinaryY" = "GB", "GLMM-native" = "NAT", "GLMM-AR1" = "AR1")[dgp]]
W2 <- dcast(A, model + N + ntime + label + method + sim_id ~ dgp, value.var = "estimate")
C_na <- corr_tab(W2[, .(model, N, ntime, label, method, sim_id, NAT, AR1)], "NAT", "AR1")[, pair := "native_vs_AR1"]
C_gn <- corr_tab(W2[, .(model, N, ntime, label, method, sim_id, GB, NAT)], "GB", "NAT")[, pair := "GenBinaryY_vs_native"]
C_ga <- corr_tab(W2[, .(model, N, ntime, label, method, sim_id, GB, AR1)], "GB", "AR1")[, pair := "GenBinaryY_vs_AR1"]
C2 <- rbind(C_na, C_gn, C_ga)
fwrite(C2, file.path(TABDIR, "V4_S2_corr_mechanisms.csv"))
S2 <- C2[n_pairs >= 100, .(cells = .N, median_rho = median(rho), q10 = quantile(rho, 0.1), q90 = quantile(rho, 0.9)),
         by = .(pair, method, m14 = model <= 4)][order(method, pair, m14)]
cat("Q5.2 between-mechanism rho (cells with >= 100 pairs):\n"); print(S2, digits = 3)
for (mt in c("gee", "glmer")) {
  a <- C2[pair == "native_vs_AR1" & model <= 4 & method == mt & n_pairs >= 100, median(rho)]
  b <- C2[pair == "GenBinaryY_vs_native" & model <= 4 & method == mt & n_pairs >= 100, median(rho)]
  bb <- C2[pair == "GenBinaryY_vs_native" & method == mt & n_pairs >= 100, median(rho)]
  cat(sprintf("Q5.2 %s: median rho native-AR1 (M1-M4) %.3f vs GenBinaryY-native (M1-M4) %.3f -> pass %s; |median GB-native| (all M) %.3f < 0.3: %s\n",
              mt, a, b, a > b, bb, abs(bb) < 0.3))
}
