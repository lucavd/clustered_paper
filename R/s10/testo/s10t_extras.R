# s10t_extras.R -- S10 "Testo": sign-flip test-inversion interval (S7d2 rows), adaptive quadrature (S7g rows),
# design checks from the replicate records. Definitions from Sections 2.2, 2.5, 3.7, 3.10 of the manuscript.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
res <- list(); put <- function(k, v, how) res[[k]] <<- data.table(item = k, value = v, how = how)
cp_lo <- function(k, n) ifelse(k > 0, qbeta(0.025, k, n - k + 1), 0); cp_hi <- function(k, n) ifelse(k < n, qbeta(0.975, k + 1, n - k), 1)
# ---- replicate records: design checks --------------------------------------------------------------------------
R <- qs_read(file.path(OUT, "x_recs.qs2"))
R[, est_design := generated & !design_degenerate]; R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
s2 <- R[ntime == 2 & generated == TRUE, .(p_single = 100 * mean((2 * N - n_obs) / N)), by = .(dgp, scen, N)]
put("t2_single_min_cell", s2[, min(p_single)], "min over (dgp, scen, N) of mean share of subjects with one observation at T = 2 (2N - n_obs)/N")
put("t2_single_max_cell", s2[, max(p_single)], "max of the same")
s2b <- R[ntime == 2 & generated == TRUE, .(p_single = 100 * mean((2 * N - n_obs) / N)), by = .(N)]
put("t2_single_min_N", s2b[, min(p_single)], "min over N (pooled over dgp, scen)"); put("t2_single_max_N", s2b[, max(p_single)], "max over N")
put("t2_single_all", R[ntime == 2 & generated == TRUE, 100 * sum(2 * N - n_obs) / sum(N)], "pooled share")
for (tt in c(2, 10)) for (g in c("GenBinaryY", "GLMM-native", "GLMM-AR1", "Copula-AR1")) {
  x <- R[dgp == g & scen == "M1" & ntime == tt & label == "mmlb" & generated == TRUE]
  put(sprintf("prev_M1_T%d_H1_%s", tt, g), x[, sum(sum_y) / sum(n_obs)], "pooled realised prevalence, M1, H1")
  x <- R[dgp == g & scen == "M1" & ntime == tt & label == "manual" & generated == TRUE]
  put(sprintf("prev_M1_T%d_H0_%s", tt, g), x[, sum(sum_y) / sum(n_obs)], "pooled realised prevalence, M1, H0") }
for (g in c("GenBinaryY", "GLMM-native", "GLMM-AR1")) put(paste0("prev_M1c_", g), R[dgp == g & scen == "M1c" & generated == TRUE, sum(sum_y) / sum(n_obs)], "pooled realised prevalence, M1c, both arms")
dg <- R[, .(deg = sum(design_degenerate)), by = .(dgp, grp = fifelse(scen %in% c("M6", "M7"), scen, fifelse(scen == "M1c", "M1c", "main")), label)]
print(dg); fwrite(dg, file.path(OUT, "degenerate_counts.csv"))
ct <- R[const == TRUE]; put("const_n", nrow(ct), "estimable-design datasets with constant outcome")
put("const_small_pct", 100 * ct[, mean(ntime <= 3 & N <= 25)], "share with T <= 3 and N <= 25")
put("const_n_all", R[generated == TRUE & (sum_y == 0 | sum_y == n_obs), .N], "constant outcome among all generated datasets, degenerate included")
put("datasets_n", nrow(R), "replicate records")
# ---- sign-flip test inversion (S7d2) ---------------------------------------------------------------------------
X <- qs_read(file.path(OUT, "x_S7d2.qs2"))$main
X <- merge(X[, !c("generated", "design_degenerate", "n_obs", "sum_y")], R[, c(K6, "est_design", "const"), with = FALSE], by = K6)
X[, np := NP[scen]]; X[, arm := fifelse(label == "mmlb", "H1", "H0")]; X[, small := ntime <= 3 & N <= 25]
cat("S7d2 theta by arm:\n"); print(X[, .(min = min(theta, na.rm = TRUE), max = max(theta, na.rm = TRUE), na = sum(is.na(theta))), by = .(arm, dgp)])
X[, ci_ok := est_design & !const & !is.na(lo) & !is.na(hi) & !(is.na(ci_error) == FALSE)]
X[, ci_cov := ci_ok & lo <= theta & theta <= hi]
X[, fin := ci_ok & is.finite(lo) & is.finite(hi)]
X[, g_ok := est_design & !const & ar1_gee_ok %in% TRUE & is.finite(ar1_est) & abs(ar1_est) < 10 & is.finite(ar1_se_fg) & ar1_se_fg > 0 & ar1_se_fg < 10]
X[, g_q := qt(0.975, N - np)]; X[, g_cov := g_ok & abs(ar1_est - theta) <= g_q * ar1_se_fg]
H1 <- X[arm == "H1"]
put("flipci_cov_h1_pooled", 100 * H1[, sum(ci_cov) / sum(ci_ok)], "coverage of the test-inversion interval, H1, pooled over datasets with an interval")
put("flipci_cov_h1_meandp", H1[, .(c = 100 * sum(ci_cov) / sum(ci_ok)), by = DPK][, mean(c, na.rm = TRUE)], "same, mean over design points")
put("gee_cov_h1_pooled", 100 * H1[, sum(g_cov) / sum(g_ok)], "coverage GEE AR(1) FG t, H1, pooled")
put("gee_cov_h1_meandp", H1[, .(c = 100 * sum(g_cov) / sum(g_ok)), by = DPK][, mean(c, na.rm = TRUE)], "same, mean over design points")
put("flipci_finite_small_pooled", 100 * H1[small == TRUE, sum(fin) / sum(ci_ok)], "both ends finite, H1, small designs, among datasets with an interval")
put("flipci_finite_small_meandp", H1[small == TRUE, .(f = 100 * sum(fin) / sum(ci_ok)), by = DPK][, mean(f, na.rm = TRUE)], "same, mean over design points")
put("flipci_finite_small_est", 100 * H1[small == TRUE, sum(fin) / sum(est_design)], "both ends finite, H1, small designs, per estimable dataset")
dpc <- H1[, .(n = sum(ci_ok), k = sum(ci_ok & !ci_cov)), by = c(DPK, "small")]
dpc[, lo_flag := n >= 20 & cp_hi(k, n) < 0.025]; dpc[, hi_flag := n >= 20 & cp_lo(k, n) > 0.075]
put("flipci_lo_h1", dpc[, sum(lo_flag)], "H1 design points with non-coverage clearly below 2.5% (exact interval, n >= 20)")
put("flipci_lo_h1_T2", dpc[ntime == 2, sum(lo_flag)], "of which T = 2")
put("flipci_hi_h1", dpc[, sum(hi_flag)], "non-coverage clearly above 7.5%")
W <- H1[fin & g_ok]; W[, ratio := (hi - lo) / (2 * g_q * ar1_se_fg)]
put("ci_width_pct_geo_pooled", 100 * (exp(mean(log(W$ratio))) - 1), "geometric mean width ratio flip / GEE AR(1) FG t - 1, H1, both finite, pooled")
put("ci_width_pct_geo_meandp", 100 * (exp(W[, .(m = mean(log(ratio))), by = DPK][, mean(m)]) - 1), "same, geometric mean over design points of per-point geometric means")
rm(X, H1, W); invisible(gc())
# ---- adaptive quadrature (S7g) ------------------------------------------------------------------------------------
Q <- qs_read(file.path(OUT, "x_S7g.qs2")); RQ <- Q$recs; RQ[, est_design := generated & !design_degenerate]; RQ[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
G <- merge(Q$rows[term %in% "time:trt" | is.na(term)], RQ[, c(K6, "est_design", "const"), with = FALSE], by = K6)
cat("S7g: nAGQ values", paste(sort(unique(G$nAGQ)), collapse = ","), "; dp", uniqueN(RQ[, ..DPK]), "; labels", paste(unique(RQ$label), collapse = ","), "\n")
print(unique(RQ[, .(dgp, scen)])); print(RQ[, .(Ns = paste(sort(unique(N)), collapse = ","), Ts = paste(sort(unique(ntime)), collapse = ",")), by = .(scen)])
G[, np := NP[scen]]
G[, use := est_design & !const & !is.na(term) & is.finite(estimate) & is.finite(se) & se > 0 & abs(estimate) < 10 & se < 10 & is.finite(p)]
G[, rt := use & 2 * pt(-abs(estimate / se), N - np) < 0.05]
A <- G[, .(n = sum(use), k = sum(rt)), by = c(DPK, "nAGQ")]
A <- merge(A, RQ[, .(nd = sum(est_design)), by = DPK], by = DPK)
A[, rate := 100 * k / n]; A[, `:=`(hi = n >= 20 & cp_lo(k, n) > 0.075, lo = n >= 20 & cp_hi(k, n) < 0.025)]
fwrite(A, file.path(OUT, "agq_dp.csv"))
rm(G, Q); invisible(gc())
fwrite(rbindlist(res), file.path(OUT, "extras.csv")); print(rbindlist(res))
