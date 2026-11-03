# Provenance (2026-11-02): read-only check run from a scratch directory (writes to ./out); see docs/steps/S9.md and docs/decisions/D-82.md.

# const_check.R -- read-only check: GEE fits on datasets with a constant outcome (sum_y == 0 or sum_y == n_obs)
# counted as usable. Recomputes Table 5 quantities (usable share, flags clearly above 7.5%, power) with and without them.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L)
REPO <- "/home/user/MMM"; dir.create("out", showWarnings = FALSE)
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
fd <- list.files(file.path(REPO, "results_v4/S7d/rows"), full.names = TRUE, pattern = "\\.rds$")
fe <- list.files(file.path(REPO, "results_v4/S7e/rows"), full.names = TRUE, pattern = "\\.rds$")
GC <- c("cor", "gee_ok", "est", "se_unc", "se_md", "se_kc", "se_fg", "se_cr2", "sim_id", "dgp", "scen", "N", "ntime", "label")
RC0 <- c("sim_id", "generated", "design_degenerate", "n_obs", "sum_y", "p_flip", "dgp", "scen", "N", "ntime", "label")
rd <- function(f) { x <- readRDS(f); g <- if (!is.null(x$rows)) x$rows else x$gee; list(recs = as.data.table(x$recs)[, ..RC0], gee = as.data.table(g)[, ..GC]) }
L <- parallel::mclapply(c(fd, fe), rd, mc.cores = 16L)
RECS <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); G <- rbindlist(lapply(L, `[[`, "gee"), fill = TRUE); rm(L); invisible(gc())
RECS[dgp == "copula" | grepl("^Cop", dgp), dgp := "Copula-AR1"]; G[dgp == "copula" | grepl("^Cop", dgp), dgp := "Copula-AR1"]
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
RECS[, const := generated & !design_degenerate & (sum_y == 0 | sum_y == n_obs)]
cat(sprintf("files %d + %d | recs %d (unique %d) | gee rows %d\n", length(fd), length(fe), nrow(RECS), uniqueN(RECS[, ..K6]), nrow(G)))
G <- merge(G, RECS[, c(K6, "design_degenerate", "const", "sum_y", "n_obs"), with = FALSE], by = K6)
G[, np := NP[scen]]
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
IND <- list()
for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) {
  if (ty == "cr2" && cn == "ar1") next
  X <- G[cor == cn, c(K6, "est", "np", "design_degenerate", "const", "gee_ok", SEC[[ty]]), with = FALSE]; setnames(X, SEC[[ty]], "se")
  X[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  X[, rej := use & 2 * pt(-abs(est / se), N - np) < 0.05]
  IND[[paste0("gee_", cn, "_", ty)]] <- X[, c(K6, "use", "rej", "const", "est", "se"), with = FALSE]
}
X <- RECS[, c(K6, "design_degenerate", "const", "p_flip"), with = FALSE]
X[, use := !design_degenerate & is.finite(p_flip)]; X[, rej := use & p_flip < 0.05]; IND[["flip"]] <- X[, c(K6, "use", "rej", "const"), with = FALSE]
ND <- RECS[, .(nd = sum(!design_degenerate), nconst = sum(const)), by = .(dgp, scen, N, ntime, label)]
cp_lo <- function(k, n) fifelse(k > 0, qbeta(0.025, k, n - k + 1), 0)
cp_hi <- function(k, n) fifelse(k < n, qbeta(0.975, k + 1, n - k), 1)
DP <- rbindlist(lapply(names(IND), function(s) IND[[s]][, .(n = sum(use), k = sum(rej), n_c = sum(use & const), k_c = sum(rej & const)),
                                                         by = .(dgp, scen, N, ntime, label)][, strategy := s]))
DP <- merge(DP, ND, by = c("dgp", "scen", "N", "ntime", "label"))
H0 <- DP[label == "manual"]; H1 <- DP[label == "mmlb"]
for (v in c("orig", "corr")) {
  nn <- if (v == "orig") H0$n else H0$n - H0$n_c; kk <- if (v == "orig") H0$k else H0$k - H0$k_c
  H0[, (paste0("pu_", v)) := 100 * nn / nd]
  H0[, (paste0("hi_", v)) := nn >= 20 & cp_lo(kk, nn) > 0.075]
  H0[, (paste0("lo_", v)) := nn >= 20 & cp_hi(kk, nn) < 0.025]
  H0[, (paste0("rate_", v)) := fifelse(nn > 0, 100 * kk / nn, NA_real_)]
}
H1[, `:=`(upow_orig = 100 * k / nd, upow_corr = 100 * (k - k_c) / nd, pu1_orig = 100 * n / nd, pu1_corr = 100 * (n - n_c) / nd)]
D <- merge(H0, H1[, .(dgp, scen, N, ntime, strategy, upow_orig, upow_corr, pu1_orig, pu1_corr)], by = c("dgp", "scen", "N", "ntime", "strategy"))
D[, small := ntime <= 3 & N <= 25]
fwrite(D, "out/const_dp.csv")
T5 <- D[, .(n_dp = .N, hi_orig = sum(hi_orig), hi_corr = sum(hi_corr), lo_orig = sum(lo_orig), lo_corr = sum(lo_corr),
            pu_orig = mean(pu_orig), pu_corr = mean(pu_corr), pu_small_orig = mean(pu_orig[small]), pu_small_corr = mean(pu_corr[small]),
            upow_orig = mean(upow_orig), upow_corr = mean(upow_corr), rate_orig = mean(rate_orig, na.rm = TRUE), rate_corr = mean(rate_corr, na.rm = TRUE),
            rate_small_orig = mean(rate_orig[small], na.rm = TRUE), rate_small_corr = mean(rate_corr[small], na.rm = TRUE)), by = strategy]
fwrite(T5, "out/const_table5.csv"); print(T5, digits = 4)
CNT <- rbindlist(lapply(names(IND), function(s) IND[[s]][, .(usable_const = sum(use & const), rej_const = sum(rej & const)), by = label][, strategy := s]))
fwrite(CNT, "out/const_counts.csv"); print(dcast(CNT, strategy ~ label, value.var = c("usable_const", "rej_const")))
RC <- RECS[, .(datasets = .N, nondeg = sum(!design_degenerate), const = sum(const), const_all0 = sum(const & sum_y == 0)), by = .(label, small = ntime <= 3 & N <= 25)]
fwrite(RC, "out/const_recs.csv"); print(RC)
EX <- IND[["gee_ar1_fg"]][use & const, .(est = summary(est), se = summary(se), absz = summary(abs(est / se)))]
print(EX); fwrite(EX, "out/const_ar1fg_est_summary.csv")
# baseline gate: original values must equal the paper's Table 5 file
P <- fread(file.path(REPO, "tables/V5_S7f_table4.csv"))
GT <- merge(T5[, .(strategy, hi_orig, pu_orig, pu_small_orig, upow_orig)], P[, .(strategy, hi_all, pu, pu_small, upow)], by = "strategy")
GT[, ok := hi_orig == hi_all & abs(pu_orig - pu) < 1e-9 & abs(pu_small_orig - pu_small) < 1e-9 & abs(upow_orig - upow) < 1e-9]
print(GT, digits = 6); fwrite(GT, "out/const_gate_baseline.csv")
