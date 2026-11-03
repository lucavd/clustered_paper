#!/usr/bin/env Rscript
# =============================================================================
# s7f_extras.R -- S7f: values for the v5 registry that are not a single cell of an existing CSV, and the file of the
# new comparisons (V5_S7f_claims.csv). Each value is computed from a file and written with its source and how.
# Inputs read-only. Outputs: tables/V5_S7f_extra.csv, tables/V5_S7f_claims.csv, tables/V5_S7f_rule_t2.csv,
# tables/V5_S7f_fig1c.csv, tables/V5_S7f_software.csv, results_v4/S7f/sessionInfo.txt.
# Run: Rscript --vanilla R/s7f/s7f_extras.R   (from ~/MMM; after R/s7f/s7f_table4.R)
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
setwd("/home/user/MMM")
OUT <- list()
put <- function(key, value, source, how) {
  stopifnot(length(value) == 1, !is.na(value))
  OUT[[key]] <<- data.frame(key = key, value = as.character(value), source = source, how = how, stringsAsFactors = FALSE)
}
num_in <- function(file, pattern) {
  L <- readLines(file, warn = FALSE); m <- regmatches(L, regexec(pattern, L)); m <- m[lengths(m) > 1]
  if (!length(m)) stop("pattern not found: ", pattern, " in ", file); as.numeric(m[[1]][2])
}
## ---- constants read from the executed code ----------------------------------------------------------------------
put("cop_rho", num_in("R/v4/rundir/copula/copula_params.R", "^COPULA_RHO\\s*<-\\s*([0-9.]+)"), "R/v4/rundir/copula/copula_params.R", "COPULA_RHO")
put("flip_B", num_in("R/s7c/s7c_tests.R", "B_FLIP <- ([0-9]+)L"), "R/s7c/s7c_tests.R", "B_FLIP (random sign vectors besides the identity)")
put("flip_enum_maxn", num_in("R/s7c/s7c_tests.R", "ENUM_FLIP_MAXN <- ([0-9]+)L"), "R/s7c/s7c_tests.R", "all 2^N signs up to this N")
put("flip_coefmax", num_in("R/s7c/s7c_tests.R", "^COEF_MAX <- ([0-9.]+)"), "R/s7c/s7c_tests.R", "COEF_MAX of the sensitivity filter")
put("guard_div", num_in("R/s7d2/s7d2_ci.R", "^CI_DIVERGE <- ([0-9e.+-]+)"), "R/s7d2/s7d2_ci.R", "CI_DIVERGE (D-57, D-59)")
put("guard_u0", num_in("R/s7d2/s7d2_ci.R", "^CI_UZERO <- ([0-9e.+-]+)"), "R/s7d2/s7d2_ci.R", "CI_UZERO (D-57, D-59)")
## ---- Li-Redden degrees of freedom (N - 2) on the 1,896 points of the three study mechanisms ------------------------
A <- fread("tables/V5_S7f_cp_S7d_dp_rates.csv"); C <- fread("tables/V5_S7f_cp_S7b_rule_dp.csv")   # exact-interval flags (s7f_cpflags.R)
est <- c("gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_ar1_md", "gee_ind_md", "gee_exch_md", "gee_ar1_kc", "gee_ind_kc", "gee_exch_kc",
         "gee_ind_cr2", "gee_exch_cr2", "gee_ar1_unc", "gee_ind_unc", "gee_exch_unc")
rt <- function(ref_) rbind(A[ref == ref_ & strategy %in% c(est, "study_mmm"), .(dgp, scen, N, ntime, strategy, hi, lo, pu)][strategy == "study_mmm", strategy := "mmm_cc"],
                          C[ref == ref_ & method %in% c("glmer", "shrink"), .(dgp, scen, N, ntime, strategy = method, hi, lo, pu)])
R2 <- rbind(rt("t")[, df := "N - p"], rt("t2")[, df := "N - 2"])
stopifnot(R2[, .N, by = .(strategy, df)][, all(N == 1896L)])
RU <- R2[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_lo_N10 = sum(lo & N == 10), pu = mean(pu)), by = .(df, strategy)]
RU[, eligible := n_hi == 0L]; RU[, chosen := eligible & pu == max(pu[eligible]), by = df]
setorder(RU, df, -eligible, -pu); fwrite(RU, "tables/V5_S7f_rule_t2.csv")
## ---- Figure 1 panel c: share of design points clearly above 7.5% by N, over the 2,376 points -------------------------
D <- fread("tables/V5_S7f_rule_dp.csv")
F1 <- D[strategy %in% c("gee_ar1_fg", "mmm_cc", "gee_exch_fg", "gee_exch_md", "glmer_lrt"), .(n_dp = .N, n_hi = sum(hi)), by = .(strategy, N)]
setorder(F1, strategy, N); fwrite(F1, "tables/V5_S7f_fig1c.csv")
## ---- power in the short-series corner (ceiling for the text) and width of the test-inversion interval ------------------
put("upow_small_ceiling", ceiling(D[strategy %in% c("gee_ar1_fg", "mmm_cc") & small == TRUE, .(u = mean(upow)), by = strategy][, max(u)]),
    "tables/V5_S7f_rule_dp.csv", "ceiling of the larger mean upow over T <= 3 and N <= 25, gee_ar1_fg and mmm_cc")
BM <- fread("tables/V4_S7d2_boot_means.csv")
put("ci_width_pct", 100 * (BM[q == "PG7 width ratio flip / GEE AR(1) FG t (geometric), H1", est] - 1), "tables/V4_S7d2_boot_means.csv", "100 (PG7 - 1)")
put("ci_width_pct_N75", 100 * (fread("tables/V4_S7d2_width_groups.csv")[group == "N 75-145", ratio] - 1), "tables/V4_S7d2_width_groups.csv", "100 (ratio - 1), N 75-145")
## ---- copula: model SE of the GLMM (random intercept only, M1-M2, N >= 85) over the spread of its estimates -----------
G <- fread("tables/V4_S7e_glmm_se_ratio.csv")[scen %in% c("M1", "M2") & N >= 85, .(ratio = mean(ratio)), by = .(scen, ntime)]
put("cop_glmm_ser_min", min(G$ratio), "tables/V4_S7e_glmm_se_ratio.csv", "min over (scen, T) of the mean ratio over N >= 85, M1-M2")
put("cop_glmm_ser_max", max(G$ratio), "tables/V4_S7e_glmm_se_ratio.csv", "max over (scen, T) of the mean ratio over N >= 85, M1-M2")
## ---- copula: GEE exchangeable + FG + t clearly above 7.5%: occasions and cluster range -------------------------------
B <- fread("tables/V5_S7f_cp_S7e_dp_rates.csv")[strategy == "gee_exch_fg" & ref == "t" & hi == TRUE]
stopifnot(nrow(B) > 0); tt <- B[, .N, by = ntime][order(-N)]
put("cop_exfg_T", tt$ntime[1], "tables/V5_S7f_cp_S7e_dp_rates.csv", "most frequent ntime of the points clearly above, gee_exch_fg t")
put("cop_exfg_nT", tt$N[1], "tables/V5_S7f_cp_S7e_dp_rates.csv", "number of those points at that ntime")
put("cop_exfg_Nmin", min(B$N), "tables/V5_S7f_cp_S7e_dp_rates.csv", "min N of those points"); put("cop_exfg_Nmax", max(B$N), "tables/V5_S7f_cp_S7e_dp_rates.csv", "max N")
## ---- marginal null under GLMM-native: largest expected excess rejection, marginal methods without the shrinkage cell ---
NS0 <- fread("tables/V4_S7b_null_summary.csv"); NS <- NS0[ref == "n" & method %in% c("mmm", "gee", "gee_fg", "gee_md", "gee_kc")]
put("null_exc_main", NS[group == "M1-M5", max(max_excess)], "tables/V4_S7b_null_summary.csv", "max max_excess, M1-M5, normal, MMM and GEEs")
put("null_exc_m67", NS[group == "M6-M7", max(max_excess)], "tables/V4_S7b_null_summary.csv", "max max_excess, M6-M7, normal, MMM and GEEs")
put("null_flips", fread("tables/V5_S7f_cp_null_flips.csv")$value, "tables/V5_S7f_cp_null_flips.csv", "flags (exact intervals) that change after subtraction, all rows")
## ---- fewclust against a computation by hand (S7d check GE4): largest difference over the three corrections -------------
GE <- fread("tables/V4_S7d_GE4.csv")
put("ge4_max", max(GE$max_md, GE$max_kc, GE$max_fg), "tables/V4_S7d_GE4.csv", "max of max_md, max_kc, max_fg over working correlations")
## ---- sensitivity: fits with warnings excluded (U1), main grid, H0, pooled replicate-weighted (as R/s7/s7_extras.R) --------
PA <- fread("tables/variant_cv/V4_perf_all.csv.gz")
U <- PA[arm == "H0" & scen != "M1c" & set == "own" & usab %in% c("U0", "U1") & n_obs > 0 & is.finite(rej)]
S13 <- U[, .(rate = 100 * sum(rej * n_obs) / sum(n_obs), n = sum(n_obs)), by = .(dgp, method, usab)]
S13 <- dcast(S13, dgp + method ~ usab, value.var = c("rate", "n"))
S13[, mcse_U1 := 100 * sqrt((rate_U1 / 100) * (1 - rate_U1 / 100) / n_U1)]
fwrite(S13, "tables/V5_S7f_u1.csv")
## ---- sizes --------------------------------------------------------------------------------------------------------------
T4 <- fread("tables/V5_S7f_table4.csv")
put("ndp_all5", T4[strategy == "gee_ar1_fg", n_dp], "tables/V5_S7f_table4.csv", "n_dp of the rule")
put("ndp_cop", T4[strategy == "gee_ar1_fg", n_copula], "tables/V5_S7f_table4.csv", "copula points")
put("ndp_three", T4[strategy == "gee_ar1_fg", n_main + n_m67], "tables/V5_S7f_table4.csv", "points of the three study mechanisms")
## ---- new comparisons: block-bootstrap intervals already computed in S7d, S7e, S7d2 (estimate, lo, hi) -------------------
cls <- function(lo, hi, m) {
  if (lo >= -m && hi <= m) return("trascurabile")
  if (lo * hi > 0 && min(abs(lo), abs(hi)) > m) return("risolta e rilevante")
  if (lo * hi > 0) return("risolta, rilevanza indeterminata")
  "non risolta"
}
CL <- list()
addc <- function(name, file, q_, margin, what) {
  X <- fread(file); r <- X[q == q_]; stopifnot(nrow(r) == 1)
  CL[[name]] <<- data.table(name = name, est = r$est, lo = r$lo, hi = r$hi, margin = margin, class = cls(r$lo, r$hi, margin), source = file, q = q_, what = what)
}
addc("wc_ind", "tables/V4_S7d_boot_means.csv", "FG ind - FG exch rate (fewclust)", 1, "type I error, FG independence minus exchangeable (pp)")
addc("wc_ar1", "tables/V4_S7d_boot_means.csv", "FG ar1 - FG exch rate (fewclust)", 1, "type I error, FG AR(1) minus exchangeable (pp)")
addc("impl_fg", "tables/V4_S7d_boot_means.csv", "fewclust FG exch - study FG rate", 1, "type I error, FG fewclust minus geesmv (pp)")
addc("impl_fg_pu", "tables/V4_S7d_boot_means.csv", "pu fewclust FG exch - study FG", 2.5, "usable share, FG fewclust minus geesmv (pp)")
addc("impl_kc", "tables/V4_S7d_boot_means.csv", "fewclust KC - study KC rate", 1, "type I error, KC fewclust minus geesmv (pp)")
addc("cop_mmm", "tables/V4_S7e_boot_means.csv", "MMM cc t rate H0, copula - study mechanisms", 1, "MMM with check, t: copula minus study mechanisms (pp)")
addc("cop_gee", "tables/V4_S7e_boot_means.csv", "GEE AR(1) FG t rate H0, copula - study mechanisms", 1, "GEE AR(1) FG t: copula minus study mechanisms (pp)")
addc("cop_glmm", "tables/V4_S7e_boot_means.csv", "GLMM Wald t rate H0, copula - study mechanisms", 1, "GLMM Wald t: copula minus study mechanisms (pp)")
addc("cop_flip", "tables/V4_S7e_boot_means.csv", "flip rate H0, copula - study mechanisms", 1, "sign-flip test: copula minus study mechanisms (pp)")
addc("ci_cov", "tables/V4_S7d2_boot_means.csv", "coverage test - coverage GEE AR(1) FG t, H1", 1, "coverage, test-inversion interval minus GEE AR(1) FG t (pp)")
# KC with a t reference at ten clusters: above Bradley's upper limit (claim C24a now refuted; interval wholly above 7.5)
K <- fread("tables/variant_cv/V4_S5_claims.csv")
for (g in c("GB", "NAT", "AR1")) { r <- K[name == paste0("C24a|", g)]
  CL[[paste0("kc10_", g)]] <- data.table(name = paste0("kc10_", g), est = r$est, lo = r$lo, hi = r$hi, margin = 7.5,
    class = fifelse(r$lo > 7.5, "confermata", fifelse(r$hi < 7.5, "smentita", "indeterminata")), source = "tables/variant_cv/V4_S5_claims.csv",
    q = paste0("C24a|", g), what = "rate of KC with t at N = 10 above 7.5% (lo > 7.5)") }
CLT <- rbindlist(CL); fwrite(CLT, "tables/V5_S7f_claims.csv"); print(CLT[, .(name, est, lo, hi, margin, class)], digits = 3)
## ---- software (Table S2) and sessionInfo --------------------------------------------------------------------------
pk <- c("MMLB", "lme4", "geepack", "geesmv", "clubSandwich", "data.table", "crew", "mirai", "Matrix")
SW <- data.table(package = c("R", pk), version = c(paste(R.version$major, R.version$minor, sep = "."), vapply(pk, function(p) as.character(packageVersion(p)), "")))
fc_desc <- read.dcf("/home/user/R/s7d_lib/fewclust/DESCRIPTION")
SW <- rbind(SW, data.table(package = c("fewclust", "fewclust_commit"), version = c(fc_desc[1, "Version"], substr(readLines("/home/user/R/s7d_lib/fewclust_commit.txt")[1], 1, 7))))
fwrite(SW, "tables/V5_S7f_software.csv")
dir.create("results_v4/S7f", showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output({ for (p in pk) suppressMessages(library(p, character.only = TRUE)); print(sessionInfo()) }), "results_v4/S7f/sessionInfo.txt")
put("git_commit", system("git rev-parse --short HEAD", intern = TRUE), "git", "HEAD at registry build")

## ---- second adversarial review (2026-10-27): values for the corrected sentences ------------------------------------
.libPaths(c("/home/user/R/s7d_lib", .libPaths())); suppressMessages(library(fewclust))
fgb <- c(formals(fewclust:::fewclust_control)$fg_bound, formals(fewclust:::gee_vcov)$fg_bound); fgb <- c(fgb, formals(geesmv::GEE.var.fg)$b)   # geesmv default; the simulation and clinical calls never pass b
stopifnot(length(unique(fgb)) == 1L, is.numeric(fgb))
put("fg_bound", fgb[1], "fewclust:::fewclust_control and fewclust:::gee_vcov (commit in tables/V5_S7f_software.csv)", "default bound b of the Fay-Graubard correction, min(b, H_ii)")
D7 <- fread("tables/V5_S7f_rule_dp.csv")
put("flip_upow_small", D7[strategy == "flip" & small == TRUE, mean(upow)], "tables/V5_S7f_rule_dp.csv", "mean usable power of the sign-flip test, T <= 3 and N <= 25")
GE <- fread("tables/V4_S7d_GE4.csv"); stopifnot(max(GE$max_md, GE$max_kc, GE$max_fg) < 1e-12)   # the text says 'below 10^-12'
T7 <- fread("rev2/T7_bias.csv")[spec == "pieno" & method == "MMM" & N == 50]
SR <- fread("rev2/T3_subsample_rerun.csv")[spec == "pieno" & method == "MMM" & N == 50 & status == "ok" & is.finite(estimate)]
stopifnot(nrow(SR) == T7$n, abs(mean(SR$estimate) - T7$mean_est) < 1e-9)
mc50 <- sd(SR$estimate) / sqrt(nrow(SR))
put("lus_bias_mmm_50_mcse", mc50, "rev2/T3_subsample_rerun.csv", "MCSE of the MMM bias at N = 50, cohort (sd of usable estimates / sqrt n)")
stopifnot(abs(T7$bias) < 2 * mc50)   # the text says 'indistinguishable from zero'

D8 <- fread("tables/V5_S7f_rule_dp.csv")[strategy == "glmer" & dgp == "GLMM-native"]
put("glmer_hi_native", D8[, sum(hi)], "tables/V5_S7f_rule_dp.csv", "GLMM Wald t clearly above 7.5% (exact interval) on GLMM-native")
put("glmer_ndp_native", nrow(D8), "tables/V5_S7f_rule_dp.csv", "GLMM-native design points")
DPF <- fread("tables/V4_S7e_dependence_profiles.csv")[mech == "copula floor=0.00 rho=0.6"]
stopifnot(DPF[lag == 2, cor] > DPF[lag == 1, cor]^2, DPF[lag == 3, cor] > DPF[lag == 1, cor]^3)   # 'decays more slowly than an AR(1) with the same first-lag value'

stopifnot(abs(as.numeric(OUT[["guard_u0"]]$value) - 1e-6) < 1e-15)   # the text writes $10^{-6}$
R <- rbindlist(OUT); fwrite(R, "tables/V5_S7f_extra.csv")
cat(sprintf("extras: %d values -> tables/V5_S7f_extra.csv\n", nrow(R))); print(R[, .(key, value)], nrows = 100)
print(RU[strategy %in% c("gee_ar1_fg", "mmm_cc", "gee_exch_fg", "gee_ind_fg"), .(df, strategy, n_hi, n_lo, n_lo_N10, pu, chosen)], digits = 4)
