#!/usr/bin/env Rscript
# =============================================================================
# s8b_mdkc_s26.R -- S8b sub-task, plan item X03: rows of Table S26 (decision rule with N - p and N - 2 degrees of freedom,
# 1,896 design points of the three study mechanisms: main grid, M1c, M6-M7) for the average of the Mancl-DeRouen and
# Kauermann-Carroll variances (SE = sqrt((SE_MD^2 + SE_KC^2) / 2), Ford and Westgate), GEE AR(1) and exchangeable.
# Why a recomputation: tables/V5_S8b_mdkc_avg.csv (R/s8b/provenance/du8_mdkc_avg.R) covers the 2,376 points of the four
# mechanisms with N - p only, while Table S26 is on 1,896 points with both N - p and N - 2. Here the same code as
# du8_mdkc_avg.R (usability rule, denominators, exact Clopper-Pearson flags, n >= 20 assessed) is applied to the
# per-replicate rows of results_v4/S7d/rows (read-only), with both references.
# Gates: (1) for md, kc, fg (AR(1), exchangeable) the per-point n, k, nd and flags equal tables/V5_S7f_cp_S7d_dp_rates.csv
#   (ref t = N - p, ref t2 = N - 2); (2) the same strategies reproduce tables/V5_S7f_rule_t2.csv (n_hi, n_lo, n_lo_N10, pu);
#   (3) the avgvar flags of the main grid + M1c and of M6-M7 equal hi_main and hi_m67 of tables/V5_S8b_mdkc_avg.csv.
# Outputs: tables/V5_S8b_sub_mdkc_s26.csv, tables/V5_S8b_sub_mdkc_keys.csv, tables/V5_S8b_sub_mdkc_checks.csv.
# Run from ~/MMM: Rscript --vanilla R/s8b/sub/s8b_mdkc_s26.R   (single process, about 2 minutes)
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM"); t0 <- Sys.time()
CHK <- list(); chk <- function(name, ok, detail = "") { CHK[[length(CHK) + 1]] <<- data.table(check = name, pass = isTRUE(ok), detail = detail)
  cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
OUT <- list()
put <- function(key, value, kind, source, how, mcse = NA_real_) {
  stopifnot(length(value) == 1, !is.na(value))
  OUT[[key]] <<- data.table(key_ = key, value = as.character(value), mcse = if (is.na(mcse)) "" else as.character(mcse), kind = kind, source = source, how = how) }
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); KP <- c("dgp", "scen", "N", "ntime")
GC <- c(K6, "cor", "gee_ok", "est", "se_md", "se_kc", "se_fg")
f <- list.files("results_v4/S7d/rows", pattern = "\\.rds$", full.names = TRUE)
R <- vector("list", length(f)); G <- vector("list", length(f))
for (i in seq_along(f)) { y <- readRDS(f[i]); R[[i]] <- y$recs[, c(K6, "design_degenerate"), with = FALSE]; G[[i]] <- y$rows[cor %in% c("ar1", "exch"), ..GC] }
RECS <- rbindlist(R); GEE <- rbindlist(G); rm(R, G); invisible(gc())
cat(sprintf("read %d chunks: %d recs, %d gee rows (%.1f min)\n", length(f), nrow(RECS), nrow(GEE), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
chk("replicate keys unique", RECS[, .(cnt = .N), by = K6][, all(cnt == 1L)])
ND <- RECS[label == "manual", .(nd = sum(!design_degenerate)), by = KP]
chk("1,896 design points under the null", nrow(ND) == 1896L, as.character(nrow(ND)))
cp_lo <- function(k, n) fifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
cp_up <- function(k, n) fifelse(k >= n, 1, qbeta(0.975, k + 1, pmax(n - k, 1e-12)))
GEE[, se_avgvar := sqrt((se_md^2 + se_kc^2) / 2)]
BASE <- RECS[label == "manual", c(K6, "design_degenerate"), with = FALSE]
one <- function(cn, ty) {
  X <- merge(BASE, GEE[label == "manual" & cor == cn, c(K6, "est", "gee_ok", paste0("se_", ty)), with = FALSE], by = K6, all.x = TRUE)
  setnames(X, paste0("se_", ty), "se")
  X[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(X, "use", what = paste0("gee_", cn, "_", ty))   # D-82
  X[, `:=`(rej_t = use & 2 * pt(-abs(est / se), N - NP[scen]) < 0.05, rej_t2 = use & 2 * pt(-abs(est / se), N - 2L) < 0.05)]
  D <- merge(X[, .(n = sum(use), k_t = sum(rej_t), k_t2 = sum(rej_t2)), by = KP], ND, by = KP)
  D[, `:=`(assess = n >= 20, pu = 100 * n / nd)]
  for (r in c("t", "t2")) { k <- D[[paste0("k_", r)]]
    D[, paste0("hi_", r) := assess & 100 * cp_lo(k, n) > 7.5]; D[, paste0("lo_", r) := assess & 100 * cp_up(k, n) < 2.5] }
  D[, strategy := paste0("gee_", cn, "_", ty)][]
}
DP <- rbindlist(lapply(c("ar1", "exch"), function(cn) rbindlist(lapply(c("avgvar", "md", "kc", "fg"), function(ty) one(cn, ty)))))
DP[, grp := fifelse(scen %in% c("M6", "M7"), "m67", "main")]
## ---- gate 1: per-point equality with the study per-point file ------------------------------------------------------------
A <- fread("tables/V5_S7f_cp_S7d_dp_rates.csv")[ref %in% c("t", "t2") & strategy %in% unique(DP$strategy)]
for (r in c("t", "t2")) {
  M <- merge(DP[, c(KP, "strategy", "n", "nd", paste0("k_", r), paste0("hi_", r), paste0("lo_", r)), with = FALSE],
             A[ref == r, c(KP, "strategy", "n", "nd", "kk", "hi", "lo"), with = FALSE], by = c(KP, "strategy"), suffixes = c("", ".ref"))
  chk(sprintf("gate 1 (%s): md/kc/fg per-point n, k, nd, flags equal V5_S7f_cp_S7d_dp_rates.csv", r),
      nrow(M) == 6L * 1896L && M[, all(n == n.ref & nd == nd.ref & get(paste0("k_", r)) == kk & get(paste0("hi_", r)) == hi & get(paste0("lo_", r)) == lo)],
      sprintf("%d point x strategy rows", nrow(M)))
}
## ---- Table S26 rows ---------------------------------------------------------------------------------------------------
S26 <- rbind(DP[, .(df = "N - p", n_dp = .N, n_hi = sum(hi_t), n_lo = sum(lo_t), n_lo_N10 = sum(lo_t & N == 10), pu = mean(pu)), by = strategy],
             DP[, .(df = "N - 2", n_dp = .N, n_hi = sum(hi_t2), n_lo = sum(lo_t2), n_lo_N10 = sum(lo_t2 & N == 10), pu = mean(pu)), by = strategy])
T2 <- fread("tables/V5_S7f_rule_t2.csv")
M2 <- merge(S26, T2, by = c("df", "strategy"), suffixes = c("", ".ref"))
chk("gate 2: md/kc/fg rows reproduce V5_S7f_rule_t2.csv (both df)",
    nrow(M2) == 12L && M2[, all(n_hi == n_hi.ref & n_lo == n_lo.ref & n_lo_N10 == n_lo_N10.ref & abs(pu - pu.ref) < 1e-9)], sprintf("%d rows", nrow(M2)))
AV <- fread("tables/V5_S8b_mdkc_avg.csv")
hm <- DP[strategy %like% "avgvar", .(hi_main = sum(hi_t & grp == "main"), hi_m67 = sum(hi_t & grp == "m67")), by = strategy]
M3 <- merge(hm, AV[, .(strategy, hi_main.ref = hi_main, hi_m67.ref = hi_m67)], by = "strategy")
chk("gate 3: avgvar flags of main+M1c and M6-M7 equal V5_S8b_mdkc_avg.csv", nrow(M3) == 2L && M3[, all(hi_main == hi_main.ref & hi_m67 == hi_m67.ref)],
    paste(M3[, sprintf("%s %d+%d", strategy, hi_main, hi_m67)], collapse = "; "))
fwrite(S26[order(df, strategy)], "tables/V5_S8b_sub_mdkc_s26.csv")
## ---- keys -------------------------------------------------------------------------------------------------------------
for (s in c("gee_ar1_avgvar", "gee_exch_avgvar")) {
  for (d in c("N - p", "N - 2")) { r <- S26[strategy == s & df == d]
    for (fl in c("n_hi", "n_lo", "n_lo_N10")) put(sprintf("mdkc|%s|%s|%s", s, d, fl), r[[fl]], "int", "tables/V5_S8b_sub_mdkc_s26.csv",
        sprintf("design points of the three study mechanisms (%d), exact 95%% interval, %s df", r$n_dp, d)) }
  put(sprintf("mdkc|%s|pu", s), S26[strategy == s & df == "N - p", pu], "sig2", "tables/V5_S8b_sub_mdkc_s26.csv", "mean usable share (%) over the 1,896 points (same for both df)")
  a <- AV[strategy == s]
  put(sprintf("mdkc|%s|hi_all", s), a$hi_all, "int", "tables/V5_S8b_mdkc_avg.csv", "design points clearly above 7.5% over the 2,376 points of the four mechanisms, N - p (Table 4 basis)")
  put(sprintf("mdkc|%s|pu_all", s), a$pu, "sig2", "tables/V5_S8b_mdkc_avg.csv", "mean usable share (%) over the 2,376 points of the four mechanisms")
  put(sprintf("mdkc|%s|hi_copula", s), a$hi_copula, "int", "tables/V5_S8b_mdkc_avg.csv", "of which Copula-AR1")
}
put("mdkc|ndp_three", nrow(ND), "int", "results_v4/S7d/rows", "design points of the three study mechanisms")
KK <- rbindlist(OUT); setnames(KK, "key_", "key"); fwrite(KK, "tables/V5_S8b_sub_mdkc_keys.csv")
CK <- rbindlist(CHK); fwrite(CK, "tables/V5_S8b_sub_mdkc_checks.csv")
print(S26[order(df, strategy)], digits = 4)
cat(sprintf("keys %d | checks failed %d | %.1f min\n", nrow(KK), sum(!CK$pass), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
stopifnot(all(CK$pass))
