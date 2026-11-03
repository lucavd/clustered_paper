
# s10_lead_ccvar_check.R -- S10, lead check of variant (ii) of the convergence check on the main grid (S6cc/cc_all.rds):
# B code distribution among replacements; k/n under H0 with t(N - p) at the three main-grid points named by the text track.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L)
X <- as.data.table(readRDS("/home/user/MMM/results_v4/S6cc/cc_all.rds"))
cat("columns:", paste(names(X), collapse = ","), "\n"); cat("rows:", nrow(X), "\n")
print(X[, .N, by = .(used)])
print(X[used == "B", .N, by = .(B_code)][order(B_code)])
print(X[used == "B", .(share_B_nonconv = mean(!(B_code %in% c(1L, 2L)))), by = .(scen)])
K <- readRDS("/home/user/MMM/results_v4/S9/const_keys.rds"); K <- as.data.table(K)
pts <- data.table(dgp = c("GLMM-native", "GenBinaryY", "GenBinaryY"), scen = "M5", N = c(95L, 45L, 95L), ntime = c(2L, 3L, 2L))
np <- 5L
for (r in seq_len(nrow(pts))) { q <- pts[r]
  Y <- X[dgp == q$dgp & scen == q$scen & N == q$N & ntime == q$ntime & label == "manual"]
  kk <- K[dgp == q$dgp & scen == q$scen & N == q$N & ntime == q$ntime & label == "manual", sim_id]
  Y <- Y[!(sim_id %in% kk)]
  ev <- function(est, se) { u <- is.finite(est) & is.finite(se) & se > 0 & abs(est) < 10 & se < 10; z <- est / se
    rej <- u & 2 * pt(-abs(z), q$N - np) < 0.05; c(k = sum(rej), n = sum(u)) }
  v1 <- ev(Y$cc_est, Y$cc_se)
  useB2 <- Y$used == "B" & Y$B_code %in% c(1L, 2L)
  v2 <- ev(ifelse(useB2, Y$B_est, Y$A_est), ifelse(useB2, Y$B_se, Y$A_se))
  lo <- function(k, n) qbeta(0.025, k, n - k + 1)
  cat(sprintf("%s %s N%d T%d H0: rows %d | (i) %d/%d lo %.3f | (ii) %d/%d lo %.3f\n", q$dgp, q$scen, q$N, q$ntime, nrow(Y),
      v1["k"], v1["n"], lo(v1["k"], v1["n"]), v2["k"], v2["n"], lo(v2["k"], v2["n"])))
}
