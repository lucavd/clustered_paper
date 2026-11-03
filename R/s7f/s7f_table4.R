# =============================================================================
# s7f_table4.R -- S7f: Table 4 of manuscript v5, the rule of Figure 1 (D-50) on the 2,376 design points of the four
# mechanisms, rebuilt from the per-point files (not from V4_S7e_rule.csv), and gate GF3(c): the rebuilt counts and means
# equal V4_S7e_rule.csv for every common strategy and, for the sign-flip test, V4_S7d2_flip_rules.csv.
# Inputs : tables/V4_S7d_dp_rates.csv (1,896 points: main grid, M1c, M6-M7; S7d run, fewclust corrections),
#          tables/V4_S7e_dp_rates.csv (480 copula points), tables/V4_S7b_rule_dp.csv (GLMM Wald, shrinkage, LRT on the
#          1,896 points), tables/V4_S7d2_flip_rules_dp.csv (sign-flip test with the D-59 rules).
# Outputs: tables/V5_S7f_table4.csv (one row per strategy), tables/V5_S7f_rule_dp.csv (per point and strategy),
#          tables/V5_S7f_GF3c.csv. Read-only on every input.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
setwd("/home/user/MMM")
keep <- c("dgp", "scen", "N", "ntime", "strategy", "hi", "lo", "pu", "upow")
A  <- fread("tables/V5_S7f_cp_S7d_dp_rates.csv")
B  <- fread("tables/V5_S7f_cp_S7e_dp_rates.csv")
C  <- fread("tables/V5_S7f_cp_S7b_rule_dp.csv")
FL <- fread("tables/V5_S7f_cp_S7d2_flip_rules_dp.csv")
est <- c("gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_ar1_md", "gee_ind_md", "gee_exch_md", "gee_ar1_kc", "gee_ind_kc",
         "gee_exch_kc", "gee_ind_cr2", "gee_exch_cr2", "gee_ar1_unc", "gee_ind_unc", "gee_exch_unc")
s1 <- A[ref == "t" & strategy %in% c(est, "study_mmm"), ..keep][strategy == "study_mmm", strategy := "mmm_cc"]
s2 <- C[(ref == "t" & method %in% c("glmer", "shrink")) | (ref == "n" & method == "glmer_lrt"),
        .(dgp, scen, N, ntime, strategy = method, hi, lo, pu, upow)]
s3 <- B[(ref == "t" & strategy %in% c(est, "mmm_cc", "glmer", "shrink")) | (ref == "n" & strategy == "glmer_lrt"), ..keep]
f0 <- FL[version == "with D-57 rules" & label == "manual", .(dgp, scen, N, ntime, hi, lo, pu)]
f1 <- FL[version == "with D-57 rules" & label == "mmlb", .(dgp, scen, N, ntime, upow = pu * rate / 100)]
fl <- merge(f0, f1, by = c("dgp", "scen", "N", "ntime"))[, strategy := "flip"]
D <- rbind(s1, s2, s3, fl, use.names = TRUE)
D[, grp := fifelse(dgp == "Copula-AR1", "copula", fifelse(scen %in% c("M6", "M7"), "m67", "main"))]
D[, small := ntime <= 3 & N <= 25]
print(D[, .N, by = .(strategy, grp)][, .(n_strat = .N, sizes = paste(unique(N), collapse = "/")), by = grp])
stopifnot(D[, .N, by = strategy][, all(N == 2376L)], !anyNA(D$hi), !anyNA(D$lo), !anyNA(D$pu))
fwrite(D, "tables/V5_S7f_rule_dp.csv")
T4 <- D[, .(n_dp = .N, hi_main = sum(hi & grp == "main"), hi_m67 = sum(hi & grp == "m67"), hi_copula = sum(hi & grp == "copula"),
            hi_all = sum(hi), lo_all = sum(lo), lo_N10 = sum(lo & N == 10), n_N10 = sum(N == 10),
            pu = mean(pu), upow = mean(upow, na.rm = TRUE), pu_small = mean(pu[small]), n_small = sum(small),
            n_main = sum(grp == "main"), n_m67 = sum(grp == "m67"), n_copula = sum(grp == "copula")), by = strategy]
T4[, has_interval := !(strategy %in% c("glmer_lrt", "flip"))]
T4[, eligible := has_interval & hi_all == 0L]
T4[, chosen := eligible & pu == max(pu[eligible])]
setorder(T4, -has_interval, -eligible, -pu)
fwrite(T4, "tables/V5_S7f_table4.csv"); print(T4, digits = 4)
# ---- GF3(c) -----------------------------------------------------------------------------------------------------------
R <- fread("tables/V4_S7e_rule.csv")[set == "all four mechanisms (2,376 points)"]
M <- merge(T4[, .(strategy, hi_all, lo_all, pu, upow)], R[, .(strategy, n_hi, n_lo, mean_pu, mean_upow)], by = "strategy")
M[, ok := abs(pu - mean_pu) < 1e-9 & abs(upow - mean_upow) < 1e-9]   # flags now use exact intervals (s7f_cpflags.R); usable share and power must not move
M[, `:=`(d_hi = hi_all - n_hi, d_lo = lo_all - n_lo)]; fwrite(M, "tables/V5_S7f_GF3c_flagdiff.csv")
FR <- fread("tables/V4_S7d2_flip_rules.csv")[version == "with D-57 rules" & set == "all 2,376 points" & label == "H0"]
fl_ok <- abs(T4[strategy == "flip", pu] - FR$mean_pu) < 1e-9
G <- data.table(check = c("strategies compared with V4_S7e_rule.csv", "strategies with a usable-share or power mismatch", "same chosen strategy",
                          "sign-flip test: mean usable equal to V4_S7d2_flip_rules.csv"),
                value = c(nrow(M), M[, sum(!ok)],
                          as.integer(identical(T4[chosen == TRUE, strategy], R[chosen == TRUE, strategy])), as.integer(fl_ok)))
G[, pass := c(value[1] >= 16, value[2] == 0, value[3] == 1, value[4] == 1)]
fwrite(G, "tables/V5_S7f_GF3c.csv"); print(G); if (any(!M$ok)) print(M[ok == FALSE])
stopifnot(all(G$pass))
