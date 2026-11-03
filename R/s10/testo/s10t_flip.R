# s10t_flip.R -- S10 "Testo": sign-flip score test with the non-computability rules of Section 2.4 (p0_guard of the S7d2 rows)
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L); setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
X <- qs_read(file.path(OUT, "x_S7d2.qs2"))$main
R <- qs_read(file.path(OUT, "x_recs.qs2")); R[, est_design := generated & !design_degenerate]; R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
X <- merge(X[, c(K6, "p0_guard", "why0_guard", "p0_new", "p0_old"), with = FALSE], R[, c(K6, "est_design", "const", "p_flip"), with = FALSE], by = K6, all.y = TRUE)
cat("p0_old == recs p_flip:", X[is.finite(p0_old) & is.finite(p_flip), mean(abs(p0_old - p_flip) < 1e-12)], "\n")
print(X[est_design & !const & !is.finite(p0_guard), .N, by = why0_guard][order(-N)][1:10])
cat("constant-outcome datasets with a finite p0_guard:", X[const & is.finite(p0_guard), .N], " with finite recs p_flip:", X[const & is.finite(p_flip), .N], "\n")
X[, use := est_design & !const & is.finite(p0_guard)]; X[, rej := use & p0_guard < 0.05]
X[, use_nc := est_design & is.finite(p0_guard)]
A <- X[, .(n = sum(use), k = sum(rej), nd = sum(est_design), n_nc = sum(use_nc)), by = DPK]
fwrite(A, file.path(OUT, "flip_dp.csv"))
