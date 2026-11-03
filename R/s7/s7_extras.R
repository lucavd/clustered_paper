#!/usr/bin/env Rscript
# S7 extras: values for the manuscript registry that are not a single cell of an existing CSV.
# Each value is computed here from a file (CSV, RDS, or source code) and written to
# tables/V4_S7_extra.csv with the file it came from and how. Read-only on all inputs.
# Run: Rscript --vanilla R/s7/s7_extras.R   (from ~/MMM)
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
setwd("/home/user/MMM")
OUT <- list()
put <- function(key, value, source, how) {
  stopifnot(length(value) == 1, !is.na(value))
  OUT[[key]] <<- data.frame(key = key, value = as.character(value), source = source, how = how, stringsAsFactors = FALSE)
}
num_in <- function(file, pattern) {            # first numeric capture of `pattern` in `file`
  L <- readLines(file, warn = FALSE); m <- regmatches(L, regexec(pattern, L)); m <- m[lengths(m) > 1]
  if (!length(m)) stop("pattern not found: ", pattern, " in ", file); as.numeric(m[[1]][2])
}
txt_in <- function(file, pattern) {
  L <- readLines(file, warn = FALSE); m <- regmatches(L, regexec(pattern, L)); m <- m[lengths(m) > 1]
  if (!length(m)) stop("pattern not found: ", pattern, " in ", file); m[[1]][2]
}
vec_el <- function(s, k) as.numeric(trimws(strsplit(s, ",")[[1]]))[k]

## grid ---------------------------------------------------------------------------------------
RB <- fread("tables/V4_S2_replicates_by_cell.csv")
put("n_levels_N", uniqueN(RB$N), "tables/V4_S2_replicates_by_cell.csv", "distinct N")
put("n_levels_T", uniqueN(RB$ntime), "tables/V4_S2_replicates_by_cell.csv", "distinct ntime")
put("ndp_grid", RB[label == "mmlb", .N], "tables/V4_S2_replicates_by_cell.csv", "rows with label mmlb")
put("n_nt_comb", uniqueN(RB[, .(N, ntime)]), "tables/V4_S2_replicates_by_cell.csv", "distinct (N, ntime)")
put("nsim", max(RB$n_gen), "tables/V4_S2_replicates_by_cell.csv", "max n_gen")
put("nrep_max_short", RB[n_gen < 1000, max(n_gen)], "tables/V4_S2_replicates_by_cell.csv", "max n_gen among cells with n_gen < 1000")
PP <- fread("tables/V4_perf_primary.csv", select = c("dgp", "scen", "N", "ntime", "arm", "method"))
put("ndp_m1c", uniqueN(PP[scen == "M1c", .(dgp, N, ntime)]), "tables/V4_perf_primary.csv", "distinct (dgp, N, ntime) in scenario M1c")
NSD <- fread("tables/V4_S6b_ns_dp.csv")
put("ndp_m67", uniqueN(NSD[, .(dgp, scen, N, ntime)]), "tables/V4_S6b_ns_dp.csv", "distinct (dgp, scen, N, ntime)")
DG <- fread("tables/V4_S2_degenerate_designs.csv")
put("n_degenerate", DG[deg == TRUE, .N], "tables/V4_S2_degenerate_designs.csv", "rows with deg == TRUE (keys per arm)")

## generating parameters ----------------------------------------------------------------------
PAR <- fread("tables/V4_dgm_parameters.csv")
pv <- function(dgp_, m, arm_, par) { v <- PAR[dgp == dgp_ & model == m & arm == arm_ & parameter == par, value]; stopifnot(length(v) == 1); v }
put("p_trt", as.numeric(sub(".*prob = ([0-9.]+).*", "\\1", pv("GenBinaryY", 1, "H1", "draw_03"))), "tables/V4_dgm_parameters.csv", "prob of rbinom draw_03 (trt), GenBinaryY M1 H1")
put("beta3_m4", vec_el(pv("GenBinaryY", 4, "H1", "beta"), 4), "tables/V4_dgm_parameters.csv", "4th marginal beta, GenBinaryY M4 H1")
put("gby_sigma_m3", vec_el(pv("GenBinaryY", 3, "H1", "sigma"), 2), "tables/V4_dgm_parameters.csv", "sigma slope, GenBinaryY M3 H1")
put("m3_sd0", vec_el(pv("GenBinaryY", 3, "H1", "sigma"), 1), "tables/V4_dgm_parameters.csv", "sigma intercept, GenBinaryY M3 H1")
put("mm_q", as.numeric(pv("GenBinaryY", 1, "H1", "q")), "tables/V4_dgm_parameters.csv", "q (quadrature points), GenBinaryY; mm() uses the same q (R/analysis_model3.R)")
Dm3 <- pv("GLMM-native", 3, "H1", "D")
Dv <- as.numeric(strsplit(gsub("[^0-9.;, -]", "", Dm3), "[;,]")[[1]])
put("glmm_d11", Dv[1], "tables/V4_dgm_parameters.csv", "D[1,1], GLMM-native M3 H1")
put("glmm_d22", Dv[4], "tables/V4_dgm_parameters.csv", "D[2,2], GLMM-native M3 H1")
put("zeger_f", vec_el(pv("GLMM-native", 1, "H1", "beta_c"), 4) / vec_el(pv("GenBinaryY", 1, "H1", "beta"), 4), "tables/V4_dgm_parameters.csv", "conditional / marginal time:trt, M1")
put("ar1_rho", as.numeric(pv("GLMM-AR1", 1, "H1", "rho")), "tables/V4_dgm_parameters.csv", "rho, GLMM-AR1 M1 H1")
put("ar1_s2", as.numeric(pv("GLMM-AR1", 1, "H1", "sigma2_ar")), "tables/V4_dgm_parameters.csv", "sigma2_ar, GLMM-AR1 M1 H1")
sr <- pv("GenBinaryY", 1, "H1", "seed_rule"); stopifnot(sr == "20260805L + 1000L * model + i")
put("seed_rule", "20260805 + 1000m + i", "tables/V4_dgm_parameters.csv", "seed_rule, rewritten with m = scenario")
## M6-M7 slopes, from the scenario definition
f6 <- "R/s6b/s6b_ns_common.R"; sl <- eval(parse(text = txt_in(f6, "slope = (c\\([^)]*\\))")))
put("m6_gby", sl[1], f6, "slope vector element 1 (GenBinaryY M6)"); put("m7_gby", sl[2], f6, "element 2 (GenBinaryY M7)")
put("m6_nat", sl[3], f6, "element 3 (GLMM-native M6 slope variance)"); put("m7_nat", sl[4], f6, "element 4 (GLMM-native M7)")
## time is coded 0..T-1 in the generators (R/generation_model3.R: seq(x) - 1), so the last occasion at T = 10 is 9
stopifnot(any(grepl("seq(x) - 1", readLines("R/generation_model3.R"), fixed = TRUE)))
put("m3_sd10", vec_el(pv("GenBinaryY", 3, "H1", "sigma"), 1) + 9 * vec_el(pv("GenBinaryY", 3, "H1", "sigma"), 2), "tables/V4_dgm_parameters.csv", "sigma0 + 9 slope (time 0..9)")
put("m6_sd10", 1 + 9 * sl[1], f6, "1 + 9 slope"); put("m7_sd10", 1 + 9 * sl[2], f6, "1 + 9 slope")
L <- fread("tables/V4_S2_lus_vs_dgm.csv"); li <- function(it) { v <- L[item == it, value]; stopifnot(length(v) == 1); v }
put("prev_m1_t2", li("dgm_marginal_prevalence_M1_H1_T2"), "tables/V4_S2_lus_vs_dgm.csv", "item dgm_marginal_prevalence_M1_H1_T2")
put("m1c_prev", li("prevalence_lus_severe"), "tables/V4_S2_lus_vs_dgm.csv", "item prevalence_lus_severe")
put("m1c_var", li("var_id_modelA"), "tables/V4_S2_lus_vs_dgm.csv", "item var_id_modelA")
TG <- fread("tables/V4_S2_truth_glmmnative.csv")
put("nat_h0_min", TG[label == "manual", min(truth_marginal_weighted)], "tables/V4_S2_truth_glmmnative.csv", "min truth_marginal_weighted, H0")
put("nat_h0_max", TG[label == "manual", max(truth_marginal_weighted)], "tables/V4_S2_truth_glmmnative.csv", "max truth_marginal_weighted, H0")

## methods: constants read from the executed code ---------------------------------------------
put("geesmv_timeout", num_in("R/analysis_smallsample.R", "timeout_secs = ([0-9]+)\\)"), "R/analysis_smallsample.R", "default timeout_secs")
put("usab_max", num_in("R/s4/s4_common.R", "abs\\(estimate\\) < ([0-9]+) & se > 0 & se < 10"), "R/s4/s4_common.R", "usability filter")
put("nlm_ok", paste(eval(parse(text = txt_in("R/s6/s6_mmm_cc_run.R", "A\\$code %in% (c\\([0-9L, ]+\\))"))), collapse = " or "), "R/s6/s6_mmm_cc_run.R", "accepted nlm codes")
put("cc_thr", num_in("R/s6/s6_mmm_cc_run.R", "sm < ([0-9e.-]+)\\)"), "R/s6/s6_mmm_cc_run.R", "slope threshold")
thr <- eval(parse(text = txt_in("R/s6b/s6b_thr_build.R", "thr <- (c\\(t4 = [0-9e.-]+, t2 = [0-9e.-]+\\))")))
put("cc_thr_lo", thr[["t4"]], "R/s6b/s6b_thr_build.R", "threshold of variant t4")
put("cc_thr_hi", thr[["t2"]], "R/s6b/s6b_thr_build.R", "threshold of variant t2")
put("boot_B", nrow(readRDS("results_v4/S4/boot_idx.rds")[[1]]), "results_v4/S4/boot_idx.rds", "rows of the bootstrap index matrix")
put("min_obs_dp", num_in("py/s6/s6_figures.py", "^MIN_OBS = ([0-9]+)"), "py/s6/s6_figures.py", "MIN_OBS")
put("nl_trunc", num_in("py/s6/s6_figures.py", "ax.set_ylim\\(0, ([0-9]+)\\); ax.set_yticks"), "py/s6/s6_figures.py", "nested-loop y limit")
put("mcse_5pct", 100 * sqrt(0.05 * 0.95 / max(RB$n_gen)), "tables/V4_S2_replicates_by_cell.csv", "100 sqrt(0.05 0.95 / n), n = max n_gen")
put("mcse_5pct_2x", 200 * sqrt(0.05 * 0.95 / max(RB$n_gen)), "tables/V4_S2_replicates_by_cell.csv", "two MCSE at 5%")

## convergence check, study scenarios M3-M5 ----------------------------------------------------
CC <- readRDS("results_v4/S6cc/cc_all.rds")
u <- CC[A_ok %in% TRUE & A_has_term %in% TRUE]      # fits of the study model that returned the term (as in R/s6b/s6b_thr_build.R)
put("cc_nfits", nrow(u), "results_v4/S6cc/cc_all.rds", "rows with A_ok and A_has_term")
put("cc_share_study", 100 * mean(u$used == "B"), "results_v4/S6cc/cc_all.rds", "share of those rows with used == B (fit without slope)")

## tables: largest MCSE for captions --------------------------------------------------------------
# S8b: the tables of the paper read the variant with Kauermann-Carroll from fewclust (D-60), so their caption MCSE does too
V <- fread("tables/variant_cv/V4_S6_values.csv")[variant == "main"]
fam <- sub("\\|.*", "", V$name)
put("t2_mcse_rate", max(V[fam %in% c("rate", "mad"), se]), "tables/variant_cv/V4_S6_values.csv", "max se, rate and mad")
put("t2_mcse_wc", max(V[fam == "wc", se]), "tables/variant_cv/V4_S6_values.csv", "max se, wc")
put("t3_mcse", max(V[fam == "puT", se]), "tables/variant_cv/V4_S6_values.csv", "max se, puT")
put("t5_mcse", max(V[fam %in% c("raw", "sa"), se]), "tables/variant_cv/V4_S6_values.csv", "max se, raw and sa")
pt <- V[fam == "puT"]; pt[, Tn := as.integer(sub(".*\\|", "", name))]
put("pu_T4_min", floor(min(pt[Tn >= 4, est])), "tables/variant_cv/V4_S6_values.csv", "floor of min puT over T >= 4")

## clinical cohort -----------------------------------------------------------------------------------
LL <- as.data.table(readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds"))
nz <- LL[, .N, by = id]
put("lus_npat", nrow(nz), "ARTICOLO 2026/LURID_bis/lurid_long.rds", "distinct id")
put("lus_npt", sum(nz$N) / 10, "ARTICOLO 2026/LURID_bis/lurid_long.rds", "rows / 10 zones")
put("lus_nobs", LL[!is.na(lus_severe), .N], "ARTICOLO 2026/LURID_bis/lurid_long.rds", "rows with lus_severe observed")
put("lus_nmulti", nz[N > 10, .N], "ARTICOLO 2026/LURID_bis/lurid_long.rds", "ids with more than 10 rows")
ninv <- LL[, .(inv = any(invasive == 1)), by = id][, sum(inv)]
put("lus_ninv", ninv, "ARTICOLO 2026/LURID_bis/lurid_long.rds", "ids with any invasive == 1")
put("lus_prev", li("prevalence_lus_severe"), "tables/V4_S2_lus_vs_dgm.csv", "item prevalence_lus_severe")
put("lus_var_id", li("var_id_modelA"), "tables/V4_S2_lus_vs_dgm.csv", "item var_id_modelA")
put("lus_var_pt", li("var_pt_modelB"), "tables/V4_S2_lus_vs_dgm.csv", "item var_pt_modelB")
put("lus_chi2", li("LRT_pt_level_B_vs_A_chi2"), "tables/V4_S2_lus_vs_dgm.csv", "item LRT_pt_level_B_vs_A_chi2")
D3 <- fread("rev2/T3_degeneracy.csv")
put("lus_nsub", D3[N == min(N), .N], "rev2/T3_degeneracy.csv", "subsamples at the smallest N")
put("lus_N10", min(D3$N), "rev2/T3_degeneracy.csv", "min N"); put("lus_N50", max(D3$N), "rev2/T3_degeneracy.csv", "max N")
Ns <- sort(unique(D3$N)); put("lus_Ns", paste0(paste(head(Ns, -1), collapse = ", "), " and ", tail(Ns, 1)), "rev2/T3_degeneracy.csv", "distinct N")
put("lus_inv_none", 100 * D3[N == min(N), mean(n_inv_pat == 0)], "rev2/T3_degeneracy.csv", "share of draws at the smallest N with no invasive patient")
put("lus_phyper", dhyper(0, ninv, nrow(nz) - ninv, min(D3$N)), "ARTICOLO 2026/LURID_bis/lurid_long.rds", "dhyper(0, invasive, others, smallest N)")


## figure of power in the short-series corner (ceiling for the text) ----------------------------------
RA <- fread("tables/V4_S6b_rule_all.csv")
put("upow_T23_ceiling", ceiling(RA[strategy %in% c("mmm|t", "gee_fg|t") & src == "tutti" & sub == "T2-3 N<=25", max(upow)]), "tables/V4_S6b_rule_all.csv", "ceiling of max upow, mmm|t and gee_fg|t, T2-3 N<=25")

## software versions (Table S2) and sessionInfo -------------------------------------------------------
pk <- c("MMLB", "lme4", "geepack", "geesmv", "data.table", "crew", "mirai", "Matrix")
SW <- data.table(package = c("R", pk), version = c(paste(R.version$major, R.version$minor, sep = "."), vapply(pk, function(p) as.character(packageVersion(p)), "")))
fwrite(SW, "tables/V4_S7_software.csv")
dir.create("results_v4/S7", showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output({ for (p in pk) suppressMessages(library(p, character.only = TRUE)); print(sessionInfo()) }), "results_v4/S7/sessionInfo.txt")

## sensitivity: fits with warnings excluded (U1), main grid, H0, pooled replicate-weighted ----------------
PA <- fread("tables/variant_cv/V4_perf_all.csv.gz")
U <- PA[arm == "H0" & scen != "M1c" & set == "own" & usab %in% c("U0", "U1") & n_obs > 0 & is.finite(rej)]
S13 <- U[, .(rate = 100 * sum(rej * n_obs) / sum(n_obs), n = sum(n_obs)), by = .(dgp, method, usab)]
S13 <- dcast(S13, dgp + method ~ usab, value.var = c("rate", "n"))
S13[, mcse_U1 := 100 * sqrt((rate_U1 / 100) * (1 - rate_U1 / 100) / n_U1)]
fwrite(S13, "tables/V4_S7_u1.csv")

## commit -----------------------------------------------------------------------------------------------
put("git_commit", system("git rev-parse --short HEAD", intern = TRUE), "git", "HEAD at registry build")

R <- rbindlist(OUT)
fwrite(R, "tables/V4_S7_extra.csv")
cat(sprintf("extras: %d values -> tables/V4_S7_extra.csv\n", nrow(R)))
cat(sprintf("check cc_nfits = %s (S6b inventory: 1685204)\n", R[key == "cc_nfits", value]))
print(R[, .(key, value)], nrows = 200)
