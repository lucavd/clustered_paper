.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L); setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L)
S0 <- qs_read(file.path(OUT, "x_study_nocheck_mmm.qs2"))
C6 <- qs_read(file.path(OUT, "x_s6cc.qs2"))[, c(K6, "A_se", "A_p", "A_est"), with = FALSE]
Z <- merge(S0[scen %in% c("M3", "M4", "M5")], C6, by = K6)
cat("se equal:", Z[(abs(se - A_se) <= 1e-8 * pmax(1, abs(A_se))) | (is.na(se) & is.na(A_se)), .N], "of", nrow(Z), "\n")
print(summary(Z[!((abs(se - A_se) <= 1e-8 * pmax(1, abs(A_se))) | (is.na(se) & is.na(A_se))), se / A_se]))
R <- qs_read(file.path(OUT, "x_recs.qs2")); R[, est_design := generated & !design_degenerate]; R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
X <- merge(R[scen %in% c("M1", "M2", "M3", "M4", "M5") & dgp != "Copula-AR1", c(K6, "est_design", "const"), with = FALSE], S0, by = K6, all.x = TRUE)
X[, np := NP[scen]]
X[, use := est_design & !const & !is.na(term) & is.finite(estimate) & is.finite(se) & se > 0 & abs(estimate) < 10 & se < 10 & (is.na(p) | is.finite(p))]
X[, `:=`(rn = use & 2 * pnorm(-abs(estimate / se)) < 0.05, rt = use & 2 * pt(-abs(estimate / se), N - np) < 0.05)]
A <- X[, .(n = sum(use), kn = sum(rn), kt = sum(rt), nd = sum(est_design), n_U0 = sum(U0 %in% TRUE & est_design)), by = DPK]
fwrite(A, file.path(OUT, "nocheck_dp.csv"))
cp_lo <- function(k, n) ifelse(k > 0, qbeta(0.025, k, n - k + 1), 0)
H <- A[label == "manual" & scen %in% c("M3", "M4", "M5")]
H[, `:=`(hin = n >= 20 & cp_lo(kn, n) > 0.075, hit = n >= 20 & cp_lo(kt, n) > 0.075)]
cat("no-check (original rows), M3-M5 H0: dp", nrow(H), " hi normal", H[, sum(hin)], " (T 8-10:", H[ntime >= 8, sum(hin)], ")  hi t", H[, sum(hit)], "\n")
