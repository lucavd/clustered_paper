# du8_mdkc_avg.R -- S8 review, track Dubbi, doubt 8 / prediction PS5.
# Average of the Mancl-DeRouen (MD) and Kauermann-Carroll (KC) corrections for GEE AR(1) and exchangeable, recomputed
# from the per-replicate rows (results_v4/S7d/rows: main grid, M1c, M6-M7; results_v4/S7e/rows: Copula-AR1), with the
# usability rule, denominators, t reference (N - p df) and exact Clopper-Pearson flags of Table 4
# (R/s7d/s7d_agg.R, R/s7e/s7e_agg.R, R/s7f/s7f_cpflags.R, R/s7f/s7f_table4.R).
#   avgvar: SE = sqrt((SE_MD^2 + SE_KC^2) / 2)   (average of the two covariance estimators; primary)
#   avgse : SE = (SE_MD + SE_KC) / 2             (average of the two standard errors; sensitivity)
# Gate: the same code applied to md, kc, fg reproduces tables/V5_S7f_table4.csv exactly, and the per-point flags of
# tables/V5_S7f_rule_dp.csv. Read-only on ~/MMM. Output in the current directory. Single R process, no forks.
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(1L)
REPO <- "/home/user/MMM"; t0 <- Sys.time()
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
GC <- c(K6, "cor", "gee_ok", "est", "se_md", "se_kc", "se_fg")
rd <- function(dir, nm) {
  f <- list.files(file.path(REPO, dir), pattern = "\\.rds$", full.names = TRUE)
  R <- vector("list", length(f)); G <- vector("list", length(f))
  for (i in seq_along(f)) { y <- readRDS(f[i]); R[[i]] <- y$recs[, c(K6, "design_degenerate"), with = FALSE]
    G[[i]] <- y[[nm]][cor %in% c("ar1", "exch"), ..GC] }
  cat(sprintf("%s: %d chunks read (%.1f min)\n", dir, length(f), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(r = rbindlist(R), g = rbindlist(G), nfile = length(f))
}
a <- rd("results_v4/S7d/rows", "rows"); b <- rd("results_v4/S7e/rows", "gee")
RECS <- rbind(a$r, b$r); GEE <- rbind(a$g, b$g); rm(a, b); invisible(gc())
cat("recs", nrow(RECS), "unique", uniqueN(RECS[, ..K6]), "gee rows", nrow(GEE), "\n")
ND <- RECS[, .(nd = sum(!design_degenerate)), by = .(dgp, scen, N, ntime, label)]
cp_lo <- function(k, n) fifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
cp_up <- function(k, n) fifelse(k >= n, 1, qbeta(0.975, k + 1, pmax(n - k, 1e-12)))
GEE[, `:=`(se_avgvar = sqrt((se_md^2 + se_kc^2) / 2), se_avgse = (se_md + se_kc) / 2)]
BASE <- RECS[, c(K6, "design_degenerate"), with = FALSE]
one <- function(cn, ty) {
  X <- merge(BASE, GEE[cor == cn, c(K6, "est", "gee_ok", paste0("se_", ty)), with = FALSE], by = K6, all.x = TRUE)
  setnames(X, paste0("se_", ty), "se")
  X[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(X, "use", what = paste0("gee_", cn, "_", ty))   # D-82
  X[, rej := use & 2 * pt(-abs(est / se), N - NP[scen]) < 0.05]
  D <- X[, .(n = sum(use), kk = sum(rej)), by = .(dgp, scen, N, ntime, label)]
  D <- merge(D, ND, by = c("dgp", "scen", "N", "ntime", "label"))
  H0 <- D[label == "manual"]; H1 <- D[label == "mmlb", .(dgp, scen, N, ntime, upow = 100 * kk / nd, pu_h1 = 100 * n / nd)]
  H0[, rate := fifelse(n > 0, 100 * kk / n, NA_real_)]
  H0[, mc := fifelse(n > 0, 100 * sqrt((rate / 100) * (1 - rate / 100) / n), NA_real_)]
  H0[, `:=`(assess = n >= 20, pu = 100 * n / nd)]
  H0[, `:=`(cp_l = 100 * cp_lo(kk, n), cp_u = 100 * cp_up(kk, n))]
  H0[, `:=`(hi = assess & cp_l > 7.5, lo = assess & cp_u < 2.5, hi_wald = assess & rate - 2 * mc > 7.5)]
  for (v in c("hi", "lo", "hi_wald")) H0[is.na(get(v)), (v) := FALSE]
  merge(H0, H1, by = c("dgp", "scen", "N", "ntime"))[, strategy := paste0("gee_", cn, "_", ty)]
}
DP <- rbindlist(lapply(c("ar1", "exch"), function(cn) rbindlist(lapply(c("avgvar", "avgse", "md", "kc", "fg"), function(ty) one(cn, ty)))))
DP[, grp := fifelse(dgp == "Copula-AR1", "copula", fifelse(scen %in% c("M6", "M7"), "m67", "main"))]
DP[, small := ntime <= 3 & N <= 25]
fwrite(DP, "du8_dp_mdkc_avg.csv")
T4 <- DP[, .(n_dp = .N, hi_main = sum(hi & grp == "main"), hi_m67 = sum(hi & grp == "m67"), hi_copula = sum(hi & grp == "copula"),
             hi_all = sum(hi), lo_all = sum(lo), lo_N10 = sum(lo & N == 10), n_N10 = sum(N == 10),
             pu = mean(pu), upow = mean(upow, na.rm = TRUE), pu_small = mean(pu[small]), n_small = sum(small),
             n_main = sum(grp == "main"), n_m67 = sum(grp == "m67"), n_copula = sum(grp == "copula"),
             hi_wald_all = sum(hi_wald), n_assessed = sum(assess), mean_rate = mean(rate, na.rm = TRUE), max_rate = max(rate[assess], na.rm = TRUE),
             mad = mean(abs(rate - 5), na.rm = TRUE), n_point_above = sum(assess & rate > 7.5)), by = strategy]
fwrite(T4, "du8_table4_mdkc_avg.csv"); print(T4, digits = 4)
# ---- gate: reproduce Table 4 rows of md, kc, fg ------------------------------------------------------------------
REF <- fread(file.path(REPO, "tables/V5_S7f_table4.csv"))
cols <- c("n_dp", "hi_main", "hi_m67", "hi_copula", "hi_all", "lo_all", "lo_N10", "n_N10", "pu", "upow", "pu_small", "n_small", "n_main", "n_m67", "n_copula")
M <- merge(T4[, c("strategy", cols), with = FALSE], REF[, c("strategy", cols), with = FALSE], by = "strategy", suffixes = c("", ".ref"))
M[, maxdiff := apply(abs(as.matrix(M[, cols, with = FALSE]) - as.matrix(M[, paste0(cols, ".ref"), with = FALSE])), 1, max)]
G <- M[, .(strategy, maxdiff, pass = maxdiff < 1e-9)]; print(G); fwrite(G, "du8_gate_table4.csv")
RD <- fread(file.path(REPO, "tables/V5_S7f_rule_dp.csv"))
P <- merge(DP[strategy %in% REF$strategy, .(dgp, scen, N, ntime, strategy, hi, lo, pu, upow)],
           RD[, .(dgp, scen, N, ntime, strategy, hi.r = hi, lo.r = lo, pu.r = pu, upow.r = upow)], by = c("dgp", "scen", "N", "ntime", "strategy"))
GP <- P[, .(n_points = .N, hi_diff = sum(hi != hi.r), lo_diff = sum(lo != lo.r), pu_maxdiff = max(abs(pu - pu.r)),
            upow_maxdiff = max(abs(upow - upow.r), na.rm = TRUE)), by = strategy]
print(GP); fwrite(GP, "du8_gate_points.csv")
# ---- the flagged points of the averages ----------------------------------------------------------------------------
HI <- DP[strategy %like% "avg" & hi == TRUE, .(strategy, dgp, scen, N, ntime, n, kk, rate, cp_l, cp_u)][order(strategy, -cp_l)]
fwrite(HI, "du8_hi_points_avg.csv"); print(HI[, .N, by = .(strategy, dgp)])
print(head(HI[strategy == "gee_ar1_avgvar"], 20)); print(head(HI[strategy == "gee_exch_avgvar"], 10))
cat(sprintf("done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
