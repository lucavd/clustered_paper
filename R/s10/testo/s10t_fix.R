# s10t_fix.R -- S10 "Testo": (1) adaptive quadrature: paired AGQ10 - Laplace difference, with and without the
# constant-outcome rule; usable constant-outcome fits under the S7g usability; (2) width of the test-inversion interval.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
res <- list(); put <- function(k, v, how) res[[length(res) + 1]] <<- data.table(item = k, value = v, how = how)
Q <- qs_read(file.path(OUT, "x_S7g.qs2")); RQ <- Q$recs; RQ[, est_design := generated & !design_degenerate]; RQ[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
G <- merge(Q$rows[nAGQ %in% c(1L, 10L)], RQ[, c(K6, "est_design", "const"), with = FALSE], by = K6); G[, np := NP[scen]]
G[, use_nc := est_design & !is.na(term) & converged %in% TRUE & is.finite(estimate) & is.finite(se) & se > 0 & abs(estimate) < 10 & se < 10 & is.finite(p)]
G[, use := use_nc & !const]; G[, rt := 2 * pt(-abs(estimate / se), N - np) < 0.05]
put("agq_const_usable_rows", G[use_nc & const, .N], "S7g fits usable under the S7g rule although the outcome is constant (nAGQ 1 and 10)")
W <- dcast(G, dgp + scen + N + ntime + label + sim_id ~ nAGQ, value.var = c("use", "use_nc", "rt"))
for (v in c("use", "use_nc")) {
  b <- W[get(paste0(v, "_1")) & get(paste0(v, "_10"))]
  D <- b[, .(d = 100 * mean(as.numeric(rt_10) - as.numeric(rt_1)), n = .N), by = .(dgp, scen, N, ntime)]
  M <- D[, .(mean_diff = mean(d)), by = dgp]
  for (i in seq_len(nrow(M))) put(sprintf("agq_paired_%s_%s", v, M$dgp[i]), M$mean_diff[i], "mean over design points of the paired difference AGQ10 - Laplace (t), replicates usable for both")
}
X <- qs_read(file.path(OUT, "x_S7d2.qs2"))$main
R <- qs_read(file.path(OUT, "x_recs.qs2")); R[, est_design := generated & !design_degenerate]; R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
X <- merge(X[label == "mmlb", !c("generated", "design_degenerate", "n_obs", "sum_y")], R[, c(K6, "est_design", "const"), with = FALSE], by = K6)
X[, np := NP[scen]]
X[, ok := est_design & !const & !is.na(lo) & !is.na(hi) & is.na(ci_error) & is.finite(lo) & is.finite(hi)]
X[, g_ok := est_design & !const & ar1_gee_ok %in% TRUE & is.finite(ar1_est) & abs(ar1_est) < 10 & is.finite(ar1_se_fg) & ar1_se_fg > 0 & ar1_se_fg < 10]
W <- X[ok & g_ok]; W[, wf := hi - lo]; W[, wg := 2 * qt(0.975, N - np) * ar1_se_fg]; W[, ratio := wf / wg]
put("width_pairs", nrow(W), "H1 datasets with both intervals finite"); put("width_zero", W[wf <= 0, .N], "of which test-inversion width 0")
V <- W[wf > 0]
put("ci_width_pct_geo_pooled", 100 * (exp(mean(log(V$ratio))) - 1), "geometric mean ratio - 1, pooled, width > 0")
put("ci_width_pct_geo_meandp", 100 * (exp(V[, .(m = mean(log(ratio))), by = DPK][, mean(m)]) - 1), "geometric mean over points of point geometric means - 1")
put("ci_width_pct_ratio_of_means", 100 * (V[, mean(wf) / mean(wg)] - 1), "ratio of mean widths - 1")
put("ci_width_pct_median", 100 * (median(V$ratio) - 1), "median ratio - 1")
put("ci_width_pct_geo_N75", 100 * (exp(mean(log(V[N >= 75]$ratio))) - 1), "geometric, N 75-145")
fwrite(rbindlist(res), file.path(OUT, "fix.csv")); print(rbindlist(res))
