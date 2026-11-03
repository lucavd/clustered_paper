# S8 track Discussione, round 2: recompute from the per-replicate rows
# (A) coverage of the interval obtained by inverting the sign-flip test (results_v4/S7d2/rows, 47,520 chunks)
# (B) exchangeable Fay-Graubard (fewclust) with t on Copula-AR1 under H0 (results_v4/S7e/rows, 9,600 chunks)
# Definitions copied from R/s7d2/s7d2_agg.R (lines 44-49) and R/s7e/s7e_agg.R (lines 51, 71-72).
# Usage: OMP_NUM_THREADS=1 Rscript --vanilla disc_r2.R ~/MMM_rev_S8/Discussione   (one process; reads ~/MMM only)
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(1L)
args <- commandArgs(TRUE); out <- if (length(args)) path.expand(args[1]) else "."
M <- path.expand("~/MMM")
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
cp <- function(k, n) { lo <- ifelse(k == 0, 0, qbeta(0.025, k, n - k + 1)); hi <- ifelse(k == n, 1, qbeta(0.975, k + 1, n - k)); list(lo, hi) }
NPc <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)

## (A) ---------------------------------------------------------------------------------------------------------------
fs <- list.files(file.path(M, "results_v4/S7d2/rows"), full.names = TRUE, pattern = "\\.rds$")
cols <- c("dgp", "scen", "N", "ntime", "label", "sim_id", "generated", "design_degenerate", "theta", "p_theta", "est", "est_why", "lo", "hi")
X <- rbindlist(lapply(fs, function(f) readRDS(f)[, ..cols]))
X[, nd := !design_degenerate]
X[, ut := nd & is.finite(p_theta)]; X[, cov_t := ut & p_theta >= 0.05]
X[, ui := nd & is.finite(est) & is.na(est_why)]; s9_drop_const(X, "ui", what = "interval of the test"); X[, cov_i := ui & lo <= theta & theta <= hi]
X[, fin := ui & is.finite(lo) & is.finite(hi)]; X[, unb := ui & !fin]; X[, unb2 := ui & !is.finite(lo) & !is.finite(hi)]
D <- X[, .(n_rows = .N, nd = sum(nd), n_ut = sum(ut), k_cov_t = sum(cov_t), n_ui = sum(ui), k_cov_i = sum(cov_i),
           n_fin = sum(fin), k_cov_fin = sum(cov_i & fin), n_unb = sum(unb), k_cov_unb = sum(cov_i & unb), n_unb_both = sum(unb2)),
       by = .(dgp, scen, N, ntime, label)]
D[, `:=`(cov_t = 100 * k_cov_t / n_ut, cov_i = 100 * k_cov_i / n_ui, cov_fin = 100 * k_cov_fin / n_fin, share_unb = 100 * n_unb / n_ui)]
f <- function(k, n) { z <- cp(n - k, n); list(n >= 20 & z[[1]] > 0.075, n >= 20 & z[[2]] < 0.025) }   # flags on NON-coverage
D[, c("hi_t", "lo_t") := f(k_cov_t, n_ut)]; D[, c("hi_i", "lo_i") := f(k_cov_i, n_ui)]; D[, c("hi_fin", "lo_fin") := f(k_cov_fin, n_fin)]
# the "finite-only" reading: unbounded intervals counted as not covering (a reader's alternative), flags on the same n_ui
D[, c("hi_i_unbnc", "lo_i_unbnc") := f(k_cov_fin, n_ui)]
P <- fread(file.path(M, "tables/V4_S7d2_dp.csv"))
D <- merge(D, P[, .(dgp, scen, N, ntime, label, f_n_ut = n_ut, f_cov_t = cov_t, f_n_ui = n_ui, f_cov_i = cov_i, f_n_fin = n_fin)],
           by = c("dgp", "scen", "N", "ntime", "label"), all.x = TRUE)
D[, match_file := n_ut == f_n_ut & n_ui == f_n_ui & n_fin == f_n_fin & abs(cov_t - f_cov_t) < 1e-9 & abs(cov_i - f_cov_i) < 1e-9]
fwrite(D, file.path(out, "disc_r2_flip_dp.csv"))
S <- D[, .(points = .N, match_file = sum(match_file, na.rm = TRUE), assessed_i = sum(n_ui >= 20),
           hi_t = sum(hi_t), lo_t = sum(lo_t), hi_i = sum(hi_i), lo_i = sum(lo_i), hi_fin = sum(hi_fin), lo_fin = sum(lo_fin),
           hi_i_unbnc = sum(hi_i_unbnc), lo_i_unbnc = sum(lo_i_unbnc), mean_cov_t = mean(cov_t, na.rm = TRUE), mean_cov_i = mean(cov_i, na.rm = TRUE),
           mean_cov_fin = mean(cov_fin, na.rm = TRUE)), by = label]
L <- D[lo_i == TRUE, .(points = .N, T2 = sum(ntime == 2), T3 = sum(ntime == 3), N_le25 = sum(N <= 25),
                       median_share_unb = median(share_unb), min_share_unb = min(share_unb), max_share_unb = max(share_unb),
                       pooled_share_unb = 100 * sum(n_unb) / sum(n_ui), pooled_cov_i = 100 * sum(k_cov_i) / sum(n_ui),
                       pooled_cov_fin = 100 * sum(k_cov_fin) / sum(n_fin), pooled_cov_t = 100 * sum(k_cov_t) / sum(n_ut),
                       still_lo_if_unb_not_cover = sum(lo_i_unbnc), hi_if_unb_not_cover = sum(hi_i_unbnc), lo_fin = sum(lo_fin), lo_t = sum(lo_t)), by = label]
fwrite(S, file.path(out, "disc_r2_flip_summary.csv")); fwrite(L, file.path(out, "disc_r2_flip_lo_points.csv"))
set.seed(20261101)
s1 <- D[label == "mmlb" & lo_i == TRUE][sample(.N, 10)]; s0 <- D[label == "mmlb" & lo_i == FALSE & n_ui >= 20][sample(.N, 10)]
fwrite(rbind(s1[, set := "flagged"], s0[, set := "not flagged"]), file.path(out, "disc_r2_flip_sample20.csv"))

## (B) ---------------------------------------------------------------------------------------------------------------
gs <- list.files(file.path(M, "results_v4/S7e/rows"), full.names = TRUE, pattern = "_manual_c[0-9]+\\.rds$")
G <- rbindlist(lapply(gs, function(f) { z <- readRDS(f)
  g <- z$gee[cor == "exch", .(sim_id, gee_ok, est, se = se_fg)]; r <- z$recs[, .(sim_id, design_degenerate, dgp, scen, N, ntime)]
  merge(r, g, by = "sim_id", all.x = TRUE) }))
G[, np := NPc[scen]]
G[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
G[, label := "manual"]; s9_drop_const(G, "use", what = "copula GEE exchangeable FG")   # D-82 (H0 files only: label manual)
G[, rej := use & 2 * pt(-abs(est / se), N - np) < 0.05]
B <- G[, .(nd = sum(!design_degenerate), n = sum(use), k = sum(rej)), by = .(dgp, scen, N, ntime)]
B[, rate := 100 * k / n]; B[, hi := n >= 20 & cp(k, n)[[1]] > 0.075]; B[, lo := n >= 20 & cp(k, n)[[2]] < 0.025]
R7 <- fread(file.path(M, "tables/V5_S7f_cp_S7e_dp_rates.csv"))[strategy == "gee_exch_fg" & ref == "t"]
B <- merge(B, R7[, .(scen, N, ntime, f_n = n, f_rate = rate, f_hi = hi)], by = c("scen", "N", "ntime"), all.x = TRUE)
fwrite(B, file.path(out, "disc_r2_copula_exch_fg_t.csv"))
cat("A:", nrow(X), "replicates,", nrow(D), "point-arms; B:", nrow(G), "replicates,", nrow(B), "points,", sum(B$hi), "clearly above\n")
print(S); print(L); print(B[hi == TRUE])
