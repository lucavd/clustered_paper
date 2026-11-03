# =============================================================================
# s11_keys.R -- S11: registry source for numbers of the text about separated datasets (D-89) and the convergence
# check with the rule D-90 that are not a single cell of an existing CSV (token family f:n11|<key>, read by
# py/s8b_build_spec.py). Every value is computed here from a file. Output: tables/V8_S11_keys.csv
# (key, value, mcse, kind, source, how). Run after R/s7f/s7f_h2h_boot.R and R/s11/s11_degen_tables.R.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L); setwd("/home/user/MMM")
OUT <- list(); put <- function(key, value, kind, source, how) { stopifnot(length(value) == 1, !is.na(value)); OUT[[length(OUT) + 1]] <<- setnames(data.table(k_ = key, value = as.character(value), mcse = "", kind = kind, source = source, how = how), "k_", "key") }
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
## ---- separated datasets (S10 inventory, D-89) -----------------------------------------------------------------------
CLS <- fread("tables/V7_S10_degen_classes.csv"); f <- "tables/V7_S10_degen_classes.csv"
a <- CLS[scope == "all"]; stopifnot(nrow(a) == 1)
put("sep|eligible", a$eligible, "int", f, "datasets with an estimable design and a non-constant outcome")
put("sep|S", a$S, "int", f, "separated: ML estimate of the mean model does not exist (linear programming, Konis 2007)")
put("sep|S_pct", a$pct_S, "sig2", f, "separated, % of eligible")
put("sep|Si", a$Si, "int", f, "separated with an infinite interaction coefficient (detectseparation)")
put("sep|S_pct_small", CLS[scope == "region" & region == "T2-3 & N<=25", pct_S], "sig2", f, "separated, % of eligible, T 2-3 and N <= 25")
put("sep|S_pct_N45", CLS[scope == "region" & region == "N>=45", pct_S], "sig2", f, "separated, % of eligible, N >= 45")
G3 <- fread("tables/V7_S10_degen_gd3.csv"); g3 <- function(x) G3[check == x, as.numeric(value)]
put("sep|gd3_n", g3("datasets in the sample"), "int", "tables/V7_S10_degen_gd3.csv", "datasets in the stratified sample of the separation cross-check")
put("sep|gd3_sep", g3("LP separated"), "int", "tables/V7_S10_degen_gd3.csv", "separated datasets in the stratified sample of 9,504")
put("sep|gd3_big", g3("glm: |b time:trt| > 10, among LP separated"), "int", "tables/V7_S10_degen_gd3.csv", "of which plain glm gives |interaction| > 10")
DV <- fread("results_v4/S10/degen/verify_flagged_divergence.csv")
put("sep|iter_n", nrow(DV), "int", "results_v4/S10/degen/verify_flagged_divergence.csv", "usable separated datasets at the two flagged points, refitted with 25 and 200 geese iterations")
put("sep|iter_out", DV[, sum(moved_out)], "int", "results_v4/S10/degen/verify_flagged_divergence.csv", "of which the estimate moves outward with 200 iterations")
CLI <- fread("tables/V7_S10_degen_clinical.csv")[method == "(draws)"]
put("sep|resp_Si", CLI[app == "respiratory", sum(Si)], "int", "tables/V7_S10_degen_clinical.csv", "respiratory subsamples (both mean models) with an infinite treatment contrast")
put("sep|resp_draws", CLI[app == "respiratory", sum(draws)], "int", "tables/V7_S10_degen_clinical.csv", "respiratory subsamples (both mean models)")
put("sep|neo_Si", CLI[app == "neonatal cohort", sum(Si)], "int", "tables/V7_S10_degen_clinical.csv", "cohort subsamples with an infinite contrast")
## ---- usable results on separated datasets, with the current rules (D-82, D-90) ------------------------------------
CL <- readRDS("results_v4/S10/degen/classes.rds")[elig == TRUE, c(K6, "S", "Si"), with = FALSE]
PR <- as.data.table(readRDS("results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds"))
P <- merge(PR[, c(K6, "strategy", "use", "rej"), with = FALSE], CL, by = K6, all.x = TRUE); P[is.na(Si), `:=`(S = FALSE, Si = FALSE)]
fp <- "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds + results_v4/S10/degen/classes.rds"
for (s in c("mmm_cc", "gee_ar1_fg")) {
  put(sprintf("sep|%s|usable", s), P[strategy == s, sum(use)], "int", fp, "usable results, both hypotheses, 2,376 design points")
  put(sprintf("sep|%s|usable_Si", s), P[strategy == s, sum(use & Si)], "int", fp, "of which on datasets with an infinite interaction")
  put(sprintf("sep|%s|usable_Si_pct", s), 100 * P[strategy == s, sum(use & Si) / sum(use)], "sig2", fp, "% of usable results on such datasets")
}
## the design points at which GEE AR(1) + FG + t is clearly above 7.5% (Table 5)
A <- rbind(fread("tables/V5_S7f_cp_S7d_dp_rates.csv")[ref == "t"], fread("tables/V5_S7f_cp_S7e_dp_rates.csv")[ref == "t"], fill = TRUE)
HI <- A[strategy == "gee_ar1_fg" & hi == TRUE][order(dgp)]
stopifnot(nrow(HI) == fread("tables/V5_S7f_table4.csv")[strategy == "gee_ar1_fg", hi_all])
for (j in seq_len(nrow(HI))) for (s in c("gee_ar1_fg", "mmm_cc")) {
  X <- P[strategy == s & label == "manual" & dgp == HI$dgp[j] & scen == HI$scen[j] & N == HI$N[j] & ntime == HI$ntime[j]]
  pre <- sprintf("sep|fgpt|%d|%s|", j, s)
  put(paste0(pre, "n"), X[, sum(use)], "int", fp, sprintf("usable, H0, %s %s N=%d T=%d", HI$dgp[j], HI$scen[j], HI$N[j], HI$ntime[j]))
  put(paste0(pre, "k"), X[, sum(rej)], "int", fp, "rejections")
  put(paste0(pre, "n_si"), X[, sum(use & Si)], "int", fp, "usable on datasets with an infinite interaction")
  put(paste0(pre, "k_si"), X[, sum(rej & Si)], "int", fp, "rejections on such datasets")
  put(paste0(pre, "n_ns"), X[, sum(use & !S)], "int", fp, "usable on datasets that are not separated")
  put(paste0(pre, "k_ns"), X[, sum(rej & !S)], "int", fp, "rejections on datasets that are not separated")
}
## ---- the convergence check under D-90 ----------------------------------------------------------------------------
SH <- fread("tables/V8_S11_cc_shares.csv"); fs <- "tables/V8_S11_cc_shares.csv"
for (st in c("main M3-M5", "M6-M7", "Copula-AR1 M3-M5")) for (q in c("boundary", "replaced", "no_result", "no_result_of_boundary")) {
  tag <- c(`main M3-M5` = "main", `M6-M7` = "m67", `Copula-AR1 M3-M5` = "cop")[[st]]
  put(sprintf("cc|%s|%s", tag, q), SH[set == st & arm == "H0" & scen == "all" & band == "all" & dgp == "all", get(q)], "sig2", fs, sprintf("%s, H0, %% of study fits", st))
  put(sprintf("cc|%s|%s|T23", tag, q), SH[set == st & arm == "H0" & scen == "all" & band == "T2-3" & dgp == "all", get(q)], "sig2", fs, sprintf("%s, H0, T 2-3", st))
}
m67 <- SH[set == "M6-M7" & arm == "H0" & band == "all" & scen != "all"]
put("cc|m67|no_result|min", m67[, min(no_result)], "sig2", fs, "M6-M7, H0, smallest no-result share over the four scenarios")
put("cc|m67|no_result|max", m67[, max(no_result)], "sig2", fs, "largest")
## ---- variants of the check (S10 variants i-iii3, independent script R/s10/testo/s10t_ccvar.R) ----------------------
CV <- fread("tables/V7_S10_testo_cc_variants.csv"); fc <- "tables/V7_S10_testo_cc_variants.csv"
for (v in c("i", "ii", "iii", "ii3", "iii3")) {
  for (q in c("hi_all", "hi_main", "hi_m67", "hi_copula", "lo_all", "pu", "pu_small", "upow")) {
    x <- CV[variant == v & table == "Table 5 (mmm_cc, t, 2376 points)" & quantity == q, as.numeric(value)]; stopifnot(length(x) == 1)
    put(sprintf("ccv|%s|%s", v, q), x, if (startsWith(q, "hi") || startsWith(q, "lo")) "int" else "sig2", fc, sprintf("variant %s, Table 5 mmm_cc", v))
  }
  put(sprintf("ccv|%s|chosen", v), CV[variant == v & quantity == "chosen strategy", value], "text", fc, sprintf("variant %s, strategy chosen by the rule", v))
}
## ---- versions a-d of the separation sensitivity, with D-90 (R/s11/s11_degen_tables.R) -----------------------------
T5V <- fread("tables/V8_S11_degen_table5_versions.csv"); ft <- "tables/V8_S11_degen_table5_versions.csv"
for (v in c("base", "a", "b", "c", "d")) {
  for (s in c("mmm_cc", "gee_ar1_fg", "flip", "glmer_lrt")) for (q in c("hi_all", "lo_all", "pu", "pu_small", "upow")) {
    x <- T5V[version == v & strategy == s, get(q)]; stopifnot(length(x) == 1)
    put(sprintf("sepv|%s|%s|%s", v, s, q), x, if (startsWith(q, "hi") || startsWith(q, "lo")) "int" else "sig2", ft, sprintf("version %s, %s", v, s))
  }
  ch <- T5V[version == v & chosen == TRUE, strategy]; put(sprintf("sepv|%s|chosen", v), if (length(ch) == 1) ch else "none", "text", ft, sprintf("version %s, strategy chosen by the rule", v))
}
## ---- Section 3.1: mean absolute departure at T = 2, the four GEEs and the marginalized model (Table S4 values) -----
VAL <- fread("tables/variant_cv/V4_S6_values.csv")[variant == "main"]; fv <- "tables/variant_cv/V4_S6_values.csv"
mg <- VAL[name %in% sprintf("madT|T2|%s", c("mmm", "gee", "gee_md", "gee_kc", "gee_fg"))]; stopifnot(nrow(mg) == 5)
put("madT2|marg|min", mg[, min(est)], "sig2", fv, "smallest mean absolute departure at T = 2 among mmm, gee, gee_md, gee_kc, gee_fg")
put("madT2|marg|max", mg[, max(est)], "sig2", fv, "largest")
put("madT2|marg|argmin", mg[which.min(est), sub("madT\\|T2\\|", "", name)], "text", fv, "method with the smallest")
put("madT2|marg|argmax", mg[which.max(est), sub("madT\\|T2\\|", "", name)], "text", fv, "method with the largest")
## ---- Table 2 and Section 2.2: conditional-to-marginal factor of GLMM-native in M1c (F19) --------------------------
EX <- fread("tables/V4_S7_extra.csv"); ev <- function(x) as.numeric(EX[key == x, value])
kz <- 16 * sqrt(3) / (15 * pi); f1 <- sqrt(1 + kz^2 * 1)
stopifnot(abs(round(f1, 2) - ev("zeger_f")) < 1e-9)   # the stored factor 1.16 is sqrt(1 + k^2 D11) with D11 = 1
put("zeger_f_m1c", round(sqrt(1 + kz^2 * ev("m1c_var")), 2), "plain", "tables/V4_S7_extra.csv", "sqrt(1 + k^2 sigma^2), k = 16 sqrt(3) / (15 pi), sigma^2 = m1c_var (GLMM-native M1c)")
## ---- supplement corrections from S10 (F09, F12, F14) -------------------------------------------------------------
CW <- fread("results_v4/S10/supplemento/s10s_ciwidth.csv"); fw <- "results_v4/S10/supplemento/s10s_ciwidth.csv (R/s10/supplemento/s10s_main.R)"
put("ciw|ds", CW[label == "mmlb" & set == "all", gm_t], "sig2", fw, "inverted sign-flip interval wider than the GEE strategy's, geometric mean over datasets where both are finite, H1 (%)")
put("ciw|ds_N75", CW[label == "mmlb" & set == "N>=75", gm_t], "sig2", fw, "same, N >= 75 (%)")
GE <- fread("tables/V4_S7d_GE4.csv")
put("ge4|both_min", GE[, min(c(both_md, both_kc, both_fg))], "int", "tables/V4_S7d_GE4.csv", "smallest number of datasets on which fewclust and the hand computation both returned a standard error, over corrections and working correlations")
TO <- fread("tables/V7_S10_supplemento_geesmv_timeouts.csv")
put("to_max_prod", TO[, max(pct_to)], "sig2", "tables/V7_S10_supplemento_geesmv_timeouts.csv", "largest share of geesmv time-outs per design point and hypothesis in the production rows (%)")
## ---- constants of the text read from the code that used them ---------------------------------------------------
cin <- function(file, pattern) { x <- regmatches(readLines(file), regexpr(pattern, readLines(file), perl = TRUE)); stopifnot(length(x) >= 1); as.numeric(sub(pattern, "\\1", x[1], perl = TRUE)) }
fv1 <- "R/s10/degen/s10_degen_verify_flagged.R"
put("sep|iter_lo", cin(fv1, 'a <- fit\\(d, "ar1", ([0-9]+)L\\)'), "int", fv1, "geese iteration limit of the study (maxit, default)")
put("sep|iter_hi", cin(fv1, 'b <- fit\\(d, "ar1", ([0-9]+)L\\)'), "int", fv1, "geese iteration limit of the refit")
put("sep|iter_shift", cin(fv1, "moved_out = abs\\(est200\\) > abs\\(est25\\) \\+ ([0-9.]+)"), "plain", fv1, "an estimate 'moved outward' when |est200| > |est25| + this value")
fv2 <- "R/s10/testo/s10t_ccvar.R"
put("nlm_code3", max(eval(parse(text = sub(".*B_code %in% ([0-9]+:[0-9]+)\\].*", "\\1", grep("B_conv3 :=", readLines(fv2), value = TRUE))))), "int", fv2, "additional nlm code accepted by variants (ii') and (iii')")
put("smv_max_i", cin("R/s7d/s7d_common.R", "SMV_MAX_I <- ([0-9]+)L"), "int", "R/s7d/s7d_common.R", "geesmv run on replicates 1 to this value of each design point and hypothesis")
## ---- sensitivity to the thresholds of the check (F13; R/s11/s11_thr_compare.R) ----------------------------------
if (file.exists("tables/V8_S11_thr_summary.csv")) { TS <- fread("tables/V8_S11_thr_summary.csv")
  for (i in seq_len(nrow(TS))) put(TS$key[i], TS$value[i], "int", "tables/V8_S11_thr_summary.csv", "threshold variants with D-82 and D-90 against the primary threshold") }
## ---- fewclust 0.4.2 (branch package, commit f05c914; engine as 0.4.1, 45c7946): coherence condition rerun after D-90 ----------------------
FC <- fread("tables/V8_S11_fewclust_coherence.csv"); fcv <- function(x) as.numeric(FC[item == x, value])
put("fc|n_cmp", fcv("n_cmp"), "int", "tables/V8_S11_fewclust_coherence.csv", "datasets regenerated from the study compared with fewclust 0.4.1 (copy of MMM-package@45c7946:fewclust/dev/checks_0.4.1/out/c3_coherence.csv; fitting code unchanged in 0.4.2, f05c914)")
put("fc|n_coh", fcv("n_coh"), "int", "tables/V8_S11_fewclust_coherence.csv", "of which the coherence condition changed the reported fit")
put("fc|version", "0.4.2", "text", "/home/user/MMM-package/fewclust/DESCRIPTION@f05c914", "fewclust version")
put("fc|commit", "f05c914", "text", "git -C /home/user/MMM-package log", "closing commit of fewclust 0.4.2 on branch package")
K <- rbindlist(OUT); stopifnot(!anyDuplicated(K$key)); fwrite(K, "tables/V8_S11_keys.csv"); cat("keys:", nrow(K), "\n")
