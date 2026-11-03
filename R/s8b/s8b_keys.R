#!/usr/bin/env Rscript
# =============================================================================
# s8b_keys.R -- S8b: values for the registry that the corrections of the final revision need (plan item R01) and that are
# not a single cell of an existing CSV. Each value is computed here from a file and written with its source and how;
# no number is typed by hand. Output: tables/V5_S8b_keys.csv (key, value, mcse, kind, source, how).
# Run from ~/MMM after R/s7f/s7f_h2h_boot.R, R/s8b/s8b_clinical_tables.R and R/s8b/s8b_pb2_collapsed.R.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L); setwd("/home/user/MMM")
OUT <- list()
put <- function(key, value, kind, source, how, mcse = NA_real_) {
  stopifnot(length(value) == 1, !is.na(value))
  OUT[[key]] <<- data.table(key_ = key, value = as.character(value), mcse = if (is.na(mcse)) "" else as.character(mcse), kind = kind, source = source, how = how)
}
## ---- P03: copula, GEE exchangeable + FG + t clearly above 7.5% (as R/s7f/s7f_extras.R) --------------------------------
B <- fread("tables/V5_S7f_cp_S7e_dp_rates.csv")[strategy == "gee_exch_fg" & ref == "t" & hi == TRUE]
tt <- B[, .N, by = ntime][order(-N)]; B4 <- B[ntime == tt$ntime[1]]; BO <- B[ntime != tt$ntime[1]]
put("cop_exfg4_Nmin", min(B4$N), "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "min N of the points at the modal ntime")
put("cop_exfg4_Nmax", max(B4$N), "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "max N of the points at the modal ntime")
put("cop_exfg_oth_n", nrow(BO), "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "points at another ntime")
stopifnot(nrow(BO) %in% 1:2)   # S9 (D-82): with the constant-outcome rule there are two such points; written in order of ntime
setorder(BO, ntime, N)
put("cop_exfg_oth_T", BO$ntime[1], "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "ntime of the first such point")
put("cop_exfg_oth_N", BO$N[1], "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "N of the first such point")
if (nrow(BO) == 2L) {
  put("cop_exfg_oth2_T", BO$ntime[2], "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "ntime of the second such point")
  put("cop_exfg_oth2_N", BO$N[2], "int", "tables/V5_S7f_cp_S7e_dp_rates.csv", "N of the second such point") }
## ---- P07: M6-M7, GEE exchangeable + FG (fewclust, Table 4) clearly above 7.5% on GenBinaryY M7 ------------------------
R <- fread("tables/V5_S7f_rule_dp.csv")[strategy == "gee_exch_fg" & scen %in% c("M6", "M7") & hi %in% TRUE]
put("m67_exfg_hi", nrow(R), "int", "tables/V5_S7f_rule_dp.csv", "points of M6-M7 clearly above 7.5%, gee_exch_fg (equals f:t4|gee_exch_fg|hi_m67)")
put("m67_exfg_hi_gbm7", R[dgp == "GenBinaryY" & scen == "M7", .N], "int", "tables/V5_S7f_rule_dp.csv", "of which GenBinaryY M7")
## ---- P04, P05: short design points and degenerate designs by grid ---------------------------------------------------
X <- readRDS("results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds")
G <- X[strategy == "gee_ar1_fg", .(nrow = .N), by = .(dgp, scen, N, ntime, label)]
put("ndp_short_main", G[label == "mmlb" & scen %in% sprintf("M%d", 1:5) & dgp != "Copula-AR1" & nrow < 1000, .N], "int",
    "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "main-grid design points (H1) with fewer than 1000 generated datasets (check: equals x:ndp_short)")
put("ndp_short_m1c", G[label == "mmlb" & scen == "M1c" & nrow < 1000, .N], "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "M1c design points (H1) with fewer than 1000 generated datasets")
M <- X[strategy == "mmm_cc" & label == "manual"]
dg <- function(sc, d) M[scen %in% sc & dgp == d, sum(design_degenerate)]
for (d in c("GenBinaryY", "GLMM-native", "GLMM-AR1")) stopifnot(dg(sprintf("M%d", 1:5), d) == dg(sprintf("M%d", 1:5), "GenBinaryY"), dg("M1c", d) == dg("M1c", "GenBinaryY"))
for (d in c("GLMM-native")) stopifnot(dg("M6", d) == dg("M6", "GenBinaryY"), dg("M7", d) == dg("M7", "GenBinaryY"))
put("deg_main", dg(sprintf("M%d", 1:5), "GenBinaryY"), "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "degenerate designs per hypothesis, main grid, each mechanism (check: x:n_degenerate)")
put("deg_m1c", dg("M1c", "GenBinaryY"), "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "degenerate designs per hypothesis, M1c")
put("deg_m6", dg("M6", "GenBinaryY"), "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "degenerate designs per hypothesis, M6")
put("deg_m7", dg("M7", "GenBinaryY"), "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "degenerate designs per hypothesis, M7")
rm(X, M, G); invisible(gc())
## ---- P33: model fits over all runs -------------------------------------------------------------------------------------
FB <- fread("tables/V5_S8b_fits_by_run.csv")
put("n_fits_total_m", round(sum(FB$total) / 1e6, 1), "plain", "tables/V5_S8b_fits_by_run.csv", "sum of column total over runs, millions")
## ---- P23: interval of the sign-flip test (counts unbounded as covering) ------------------------------------------------
FI <- fread("tables/V5_S8b_flip_interval.csv")
put("flipci_cov_h1", FI[label == "mmlb", mean_cov_i], "sig2", "tables/V5_S8b_flip_interval.csv", "mean coverage of the inverted interval over design points, H1")
put("flipci_lo_h1", FI[label == "mmlb", lo_i], "int", "tables/V5_S8b_flip_interval.csv", "design points with non-coverage clearly below 2.5%, H1")
put("flipci_lo_h0", FI[label == "manual", lo_i], "int", "tables/V5_S8b_flip_interval.csv", "design points with non-coverage clearly below 2.5%, H0")
put("flipci_points", FI[label == "mmlb", points], "int", "tables/V5_S8b_flip_interval.csv", "design points")
## ---- P29: marginal true value of GLMM-native under H1 ------------------------------------------------------------------
TR <- fread("tables/V4_S4_truth_by_cell.csv")[dgp == "GLMM-native" & label == "mmlb"]
put("nat_h1_min", min(TR[scen %in% sprintf("M%d", 1:5), truth_marginal]), "sig2", "tables/V4_S4_truth_by_cell.csv", "min over scenarios M1-M5 and T")
put("nat_h1_max", max(TR[scen %in% sprintf("M%d", 1:5), truth_marginal]), "sig2", "tables/V4_S4_truth_by_cell.csv", "max over scenarios M1-M5 and T")
put("nat_h1_m1c_min", min(TR[scen == "M1c", truth_marginal]), "sig2", "tables/V4_S4_truth_by_cell.csv", "min over T, M1c")
put("nat_h1_m1c_max", max(TR[scen == "M1c", truth_marginal]), "sig2", "tables/V4_S4_truth_by_cell.csv", "max over T, M1c")
## ---- X03: average of the Mancl-DeRouen and Kauermann-Carroll variances ---------------------------------------------------
MK <- fread("tables/V5_S8b_mdkc_avg.csv")
put("mdkc_ar1_hi", MK[strategy == "gee_ar1_avgvar", hi_all], "int", "tables/V5_S8b_mdkc_avg.csv", "points clearly above 7.5%, AR(1)")
put("mdkc_exch_hi", MK[strategy == "gee_exch_avgvar", hi_all], "int", "tables/V5_S8b_mdkc_avg.csv", "points clearly above 7.5%, exchangeable")
## ---- PB2: size-adjusted power of the mixed model from collapsed standard errors ---------------------------------------
P2 <- fread("tables/V5_S8b_pb2_summary.csv")
put("pb2_mean", P2[quantity == "mean", est], "sig2", "tables/V5_S8b_pb2_summary.csv", "mean over cells of the share of usable H1 datasets rejected with |z| > 37", mcse = P2[quantity == "mean", mcse])
put("pb2_cells", P2[quantity == "mean", cells], "int", "tables/V5_S8b_pb2_summary.csv", "main-grid cells with a size-adjusted power of the mixed model")
put("pb2_cells_q37", P2[quantity == "mean", cells_q95_above_37], "int", "tables/V5_S8b_pb2_summary.csv", "cells whose size-adjusted critical value exceeds 37")
put("pb2_mean_pow", P2[quantity == "mean", mean_pow_sa], "sig2", "tables/V5_S8b_pb2_summary.csv", "mean over the same cells of the size-adjusted power of the mixed model")
## ---- PB5 / D-70: model-based and robust SE ratios of the clinical applications ------------------------------------------
P5 <- fread("tables/V5_S8b_pb5.csv")
put("pb5_neo", P5[application == "neonatal cohort", est], "sig2", "tables/V5_S8b_pb5.csv", "mean over N of ratio MMM - ratio GLMM", mcse = P5[application == "neonatal cohort", mcse])
put("pb5_resp", P5[application == "respiratory", est], "sig2", "tables/V5_S8b_pb5.csv", "mean over N of ratio MMM - ratio GLMM", mcse = P5[application == "respiratory", mcse])
RB <- fread("tables/V5_S8b_mmm_robust_seratio.csv")[spec == "pieno"]
put("rob_neo_min", min(RB[application == "neonatal cohort", ratio_fpc]), "sig2", "tables/V5_S8b_mmm_robust_seratio.csv", "robust-SE ratio of the MMM, cohort, min over N")
put("rob_neo_max", max(RB[application == "neonatal cohort", ratio_fpc]), "sig2", "tables/V5_S8b_mmm_robust_seratio.csv", "robust-SE ratio of the MMM, cohort, max over N")
put("rob_resp_min", min(RB[application == "respiratory", ratio_fpc]), "sig2", "tables/V5_S8b_mmm_robust_seratio.csv", "robust-SE ratio of the MMM, trial, min over N")
put("rob_resp_max", max(RB[application == "respiratory", ratio_fpc]), "sig2", "tables/V5_S8b_mmm_robust_seratio.csv", "robust-SE ratio of the MMM, trial, max over N")
SR <- fread("tables/V5_S7f_seratio_fpc.csv")[spec == "pieno" & method %in% c("MMM", "GLMM") & application == "respiratory"]
put("resp_mmm_min", min(SR[method == "MMM", ratio_fpc]), "sig2", "tables/V5_S7f_seratio_fpc.csv", "model-based SE ratio of the MMM, trial, min over N")
put("resp_mmm_max", max(SR[method == "MMM", ratio_fpc]), "sig2", "tables/V5_S7f_seratio_fpc.csv", "model-based SE ratio of the MMM, trial, max over N")
put("resp_glmm_min", min(SR[method == "GLMM", ratio_fpc]), "sig2", "tables/V5_S7f_seratio_fpc.csv", "SE ratio of the GLMM, trial, min over N")
put("resp_glmm_max", max(SR[method == "GLMM", ratio_fpc]), "sig2", "tables/V5_S7f_seratio_fpc.csv", "SE ratio of the GLMM, trial, max over N")
REFL <- fread("tables/V5_S8b_clin_mmm_ref.csv")
put("lus_mmm_se_rob", REFL[app == "neonatal cohort" & spec == "pieno", se_rob], "sig2", "tables/V5_S8b_clin_mmm_ref.csv", "robust SE of the MMM, full cohort")
put("resp_mmm_se_rob", REFL[app == "respiratory" & spec == "pieno", se_rob], "sig2", "tables/V5_S8b_clin_mmm_ref.csv", "robust SE of the MMM, full trial")
## ---- x:lus_lateral with the latent variable (D-70): lateral vs anterior contrast of the MMM on the full cohort ----------
LL <- fread("tables/V5_S8b_lus_lateral.csv")
put("lus_lateral_lv", LL$estimate, "sig2", "tables/V5_S8b_lus_lateral.csv", "regionlateral of the MMM with lv.formula = ~1, q = 20, full cohort")

## ---- P04: the replicates that could not be generated all have a degenerate design --------------------------------------
X <- readRDS("results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds"); K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
NG <- X[strategy == "mmm_cc"][!X[strategy == "gee_ar1_fg"], on = K6]
put("nongen_n", nrow(NG), "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "replicates with a design record and no GEE row (not generated), both hypotheses")
put("nongen_nondeg", NG[design_degenerate == FALSE, .N], "int", "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds", "of which with an estimable design (must be 0)")
stopifnot(NG[design_degenerate == FALSE, .N] == 0L); rm(X, NG); invisible(gc())
## ---- P23: occasions of the design points where the interval of the sign-flip test over-covers ---------------------------
FL <- fread("tables/V5_S8b_flip_lo_points.csv")
put("flipci_lo_h1_T2", FL[label == "mmlb", T2], "int", "tables/V5_S8b_flip_lo_points.csv", "of the H1 points with non-coverage clearly below 2.5%, those with T = 2")
stopifnot(FL[label == "mmlb", points] == FI[label == "mmlb", lo_i])
## ---- 3.11: bias of the marginalized model at 50 patients, in Monte Carlo standard errors ---------------------------------
BL <- merge(fread("tables/V5_S7f_lus_bias.csv")[spec == "pieno" & method == "MMM"], fread("tables/V5_S7f_lus_seratio.csv")[spec == "pieno" & method == "MMM", .(N, sd_est)], by = "N")
BL[, mcse := sd_est / sqrt(n)][, z := bias / mcse]
put("lus_bias_mmm_50_z", abs(BL[N == 50, z]), "sig2", "tables/V5_S7f_lus_bias.csv", "|bias| / MCSE of the MMM at N = 50 (MCSE = SD of the estimates / sqrt(n))")
put("lus_bias_mmm_10_z", abs(BL[N == 10, z]), "sig2", "tables/V5_S7f_lus_bias.csv", "|bias| / MCSE of the MMM at N = 10")

## ---- S04, S05: shares that round to 100 at two significant figures, printed to one decimal ---------------------------
F7 <- fread("tables/V4_S7d2_flip_rules.csv")[label == "H0" & version == "as in S7d" & set == "small designs (T 2-3, N <= 25)"]
stopifnot(nrow(F7) == 1L)
put("fl7_small_pu_1d", round(F7$mean_pu, 1), "plain", "tables/V4_S7d2_flip_rules.csv", "mean usable share of the sign-flip test without the rules, small designs, H0, one decimal")
AG <- fread("tables/V5_S7f_kc_clinical_agree.csv")[app == "respiratory"]
put("ag_resp_md_1d", round(AG[corr == "md", share_within_1pct], 1), "plain", "tables/V5_S7f_kc_clinical_agree.csv", "trial, MD fewclust within 1% of geesmv, one decimal")
put("ag_resp_fg_1d", round(AG[corr == "fg", share_within_1pct], 1), "plain", "tables/V5_S7f_kc_clinical_agree.csv", "trial, FG fewclust within 1% of geesmv, one decimal")
K <- rbindlist(OUT); setnames(K, "key_", "key"); fwrite(K, "tables/V5_S8b_keys.csv"); print(K[, .(key, value, mcse)])
