
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R"); setDTthreads(2L)
# Script 4: Tables 3, 4, 5 from the per-design-point counts of scripts 1-3 (definitions of Sections 2.5-2.6).
S1 <- readRDS(file.path(OUT, "s1_gee_dp.rds")); S2 <- readRDS(file.path(OUT, "s2_study_dp.rds"))
NS <- readRDS(file.path(OUT, "s3_ns.rds")); CO <- readRDS(file.path(OUT, "s3_cop.rds")); FLIP <- readRDS(file.path(OUT, "s3_flip.rds"))
REC <- readRDS(file.path(OUT, "s1_recs.rds"))
cols <- c(KD, "strategy", "n", "kn", "kt", "nd", "nconst")
ALL <- rbind(S1[, ..cols], S2[strategy != "gee_exch_kc", ..cols], NS$dp[, ..cols], CO$dp[, ..cols], FLIP$dp[, ..cols])
ALL <- ALL[, lapply(.SD, sum), by = c(KD, "strategy"), .SDcols = c("n", "kn", "kt", "nd", "nconst")]   # chunks of 100 -> design points
ALL[, grp := fifelse(dgp == "Copula-AR1", "copula", fifelse(scen %in% c("M6", "M7"), "m67", fifelse(scen == "M1c", "m1c", "main")))]
saveRDS(ALL, file.path(OUT, "s4_all_dp.rds"))
# ---- consistency of denominators across sources ---------------------------------------------------------------
cons <- ALL[, .(n_sources = uniqueN(paste(nd, nconst))), by = KD][n_sources > 1L]
cat("design points x label with inconsistent nd/nconst across sources:", nrow(cons), "\n")
# ---- Section 2.2 / 2.5 counts (degenerate designs, constant outcomes) ---------------------------------------
D0 <- unique(ALL[, c(KD, "grp", "nd", "nconst"), with = FALSE])
FR <- FLIP$recs[, .(nrec = .N, estim = sum(estim), const = sum(const)), by = KD]
FR <- merge(FR, unique(ALL[, c(KD, "grp"), with = FALSE]), by = KD)
deg <- FR[, .(datasets = 1000L * .N, degenerate = 1000L * .N - sum(estim), const = sum(const)), by = .(grp, dgp, scen = fifelse(grp == "main", "M1-M5", scen), label)]
fwrite(deg, file.path(OUT, "t_degenerate_counts.csv"))
cst <- FLIP$recs[const == TRUE]; cat("constant datasets", nrow(cst), "of", nrow(FLIP$recs) + 0, "; share with T<=3 & N<=25:", round(100 * mean(cst$ntime <= 3 & cst$N <= 25), 2), "\n")
ck <- readRDS(file.path(REPO, "results_v4/S9/const_keys.rds")); setDT(ck); ck[, scen := norm_scen(scen)]
cat("const keys project", nrow(ck), "; own not in project", nrow(fsetdiff(cst[, ..K6], ck[, ..K6])), "; project not in own", nrow(fsetdiff(ck[, ..K6], cst[, ..K6])), "\n")
# ---- Table 3 ----------------------------------------------------------------------------------------------------
T3M <- c(gee_fg = "gee_fg", gee_md = "gee_md", mmm = "mmm_cc", glmer = "glmer", shrink = "shrink", gee = "gee", glmer_lrt = "glmer_lrt", gee_kc = "gee_exch_kc",
         mmm_f = "mmm_cc_f", mmm_row = "mmm_projrow", gee_fg_fc = "gee_exch_fg", gee_md_fc = "gee_exch_md")
MG <- c("GenBinaryY" = "GB", "GLMM-native" = "NAT", "GLMM-AR1" = "AR1")
X3 <- ALL[grp == "main" & label == "manual" & dgp %in% names(MG)]
t3 <- rbindlist(lapply(names(T3M), function(m) {
  src <- X3[strategy == T3M[[m]]]; src[, kk := kn]
  src[, .(n = sum(n), k = sum(kk)), by = .(dgp, N)][, .(method = m, rate = 100 * sum(k) / sum(n), worst = max(abs(100 * k / n - 5)), mad = mean(abs(100 * k / n - 5)), nN = .N, n = sum(n)), by = dgp] }))
t3[, mech := MG[dgp]]; fwrite(t3, file.path(OUT, "t3_recomputed.csv")); print(t3[, .(method, mech, rate, worst, mad, nN, n)], digits = 4)
# block bootstrap of Table 3 (replicate index i within each scenario, joint across mechanisms, N, T and methods)
BR <- readRDS(file.path(OUT, "s2_boot_rep.rds")); set.seed(20261015); B <- 2000L
ids <- sort(unique(BR$sim_id)); cat("bootstrap: sim ids", length(ids), "\n")
W <- lapply(setNames(nm = c("M1", "M2", "M3", "M4", "M5")), function(m) t(rmultinom(B, length(ids), rep(1, length(ids)))))
bt <- list()
for (m in names(T3M)[1:8]) for (g in names(MG)) {
  Z <- BR[strategy == T3M[[m]] & dgp == g]; U <- 0; K <- 0
  for (sc in names(W)) { z <- Z[scen == sc]
    um <- as.matrix(dcast(z, sim_id ~ N, value.var = "u", fill = 0)[match(ids, sim_id)][, -1]); um[is.na(um)] <- 0
    km <- as.matrix(dcast(z, sim_id ~ N, value.var = "r", fill = 0)[match(ids, sim_id)][, -1]); km[is.na(km)] <- 0
    U <- U + W[[sc]] %*% um; K <- K + W[[sc]] %*% km }
  rN <- 100 * K / U; dev <- abs(rN - 5)
  bt[[paste(m, g)]] <- data.table(method = m, dgp = g, se_rate = sd(100 * rowSums(K) / rowSums(U)), se_worst = sd(apply(dev, 1, max)), se_mad = sd(rowMeans(dev)),
                                  rate0 = 100 * sum(Z$r) / sum(Z$u)) }
bt <- rbindlist(bt); fwrite(bt, file.path(OUT, "t3_boot_mcse.csv")); print(bt[, .(max_se_rate = max(se_rate), max_se_mad = max(se_mad), max_se_worst = max(se_worst))])
# ---- Table 4 ----------------------------------------------------------------------------------------------------
T4M <- c(glmer = "glmer", glmer_lrt = "glmer_lrt", mmm = "mmm_cc", gee_fg = "gee_fg", gee = "gee", mmm_f = "mmm_cc_f", mmm_row = "mmm_projrow", gee_fg_fc = "gee_exch_fg", gee_unc_fc = "gee_exch_unc")
X4 <- ALL[strategy %in% T4M & grp %in% c("main", "m1c", "copula")]
X4[, re := fifelse(scen %in% c("M1", "M2", "M1c"), "RI", "RS")]
v4 <- list(main_both = X4[grp == "main"], main_H0 = X4[grp == "main" & label == "manual"], main_H1 = X4[grp == "main" & label == "mmlb"],
           mainM1c_both = X4[grp %in% c("main", "m1c")], mainCop_both = X4[grp %in% c("main", "copula")])
t4 <- rbindlist(lapply(names(v4), function(v) v4[[v]][, .(pu = 100 * sum(n) / sum(nd), nd = sum(nd), ncells = .N), by = .(re, strategy, ntime)][, variant := v]))
t4[, se_binom := 100 * sqrt(pu / 100 * (1 - pu / 100) / nd)]
fwrite(t4, file.path(OUT, "t4_recomputed.csv")); print(t4[variant == "main_both" & strategy %in% c("glmer", "mmm_cc"), .(re, strategy, ntime, pu, se_binom)], digits = 4)
print(t4[, .(max_se = max(se_binom)), by = variant])
# ---- Table 5 ----------------------------------------------------------------------------------------------------
T5S <- c("mmm_cc", "gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_ar1_md", "gee_exch_md", "gee_ar1_kc", "gee_exch_kc", "gee_exch_cr2", "gee_exch_unc",
         "shrink", "glmer", "glmer_lrt", "flip", "mmm_cc_f", "mmm_projrow", "gee_ind_md", "gee_ind_kc", "gee_ind_cr2", "gee_ar1_unc", "gee_ind_unc")
# GEE strategies come from the fewclust rows (S7d, S7e): on the main grid drop the geesmv study rows of the same name
X5 <- ALL[strategy %in% T5S]
stopifnot(X5[, .(cnt = .N), by = c(KD, "strategy")][, all(cnt == 1L)], X5[, all(n <= nd)])
H0 <- X5[label == "manual"]; H1 <- X5[label == "mmlb", c(KD[1:4], "strategy", "n", "kt", "nd"), with = FALSE]; setnames(H1, c("n", "kt", "nd"), c("n1", "k1", "nd1"))
P5 <- merge(H0, H1, by = c(KD[1:4], "strategy"), all.x = TRUE)
P5[, k := kt]   # Wald: t with N - p; LRT and sign-flip: their own p-value (kn == kt)
P5[, `:=`(lcl = fifelse(k > 0, qbeta(0.025, k, n - k + 1), 0), ucl = fifelse(k < n, qbeta(0.975, k + 1, n - k), 1))]
P5[, assess := n >= 20L]
P5[, `:=`(hi = assess & lcl > 0.075, lo = assess & ucl < 0.025, hi_noA = n > 0 & lcl > 0.075, lo_noA = n > 0 & ucl < 0.025)]
P5[, `:=`(pu0 = 100 * n / nd, pu1 = 100 * n1 / nd1, upow = 100 * k1 / nd1)]
P5[, pub := (pu0 + pu1) / 2]; P5[, small := ntime <= 3 & N <= 25]
P5[, grp5 := fifelse(grp %in% c("main", "m1c"), "main", grp)]
fwrite(P5, file.path(OUT, "t5_points.csv"))
t5 <- P5[, .(n_dp = .N, hi_main = sum(hi & grp5 == "main"), hi_m67 = sum(hi & grp5 == "m67"), hi_copula = sum(hi & grp5 == "copula"), lo_all = sum(lo), lo_N10 = sum(lo & N == 10),
             n_N10 = sum(N == 10), pu = mean(pu0), pu_both = mean(pub), upow = mean(upow), pu_small = mean(pu0[small]), pu_small_both = mean(pub[small]), n_small = sum(small),
             hi_noA = sum(hi_noA), lo_noA = sum(lo_noA), n_unassessed = sum(!assess),
             mc_pu = mean(sqrt(pu0 * (100 - pu0) / nd)), mc_upow = mean(sqrt(upow * (100 - upow) / nd1)), mc_small = mean(sqrt(pu0 * (100 - pu0) / nd)[small])), by = strategy]
t5[, has_interval := !(strategy %in% c("glmer_lrt", "flip"))]
t5[, eligible := has_interval & (hi_main + hi_m67 + hi_copula) == 0L & !(strategy %in% c("mmm_cc_f", "mmm_projrow"))]
t5[, chosen := eligible & pu == max(pu[eligible])]
fwrite(t5, file.path(OUT, "t5_recomputed.csv")); print(t5, digits = 4)
cat("done\n")
