#!/usr/bin/env Rscript
# =============================================================================
# s8b_table_s37.R -- S8b sub-task, plan item X02: Table S37 of the supplement (Section S11, "What the decision rule
# guarantees"): how stable the count of design points clearly above 7.5% is, for the strategies of Table 4 that are in
# tables/V5_S8b_rule_stability.csv (written by R/s8b/provenance/lit_pow_rule_stability_v2.R, read-only here):
#   (1) block bootstrap over replicates (B resamples of the replicate index within model, common to all points of a
#       model): share of resamples with zero flags, mean and 90th percentile of the number of flags;
#   (2) binomial plug-in: points independent, true rate = observed rate; expected flags and P(zero flags);
#   (3) the same plug-in with rates shrunk by empirical Bayes within strategy x mechanism (method of moments).
# Also checks the ranking statement of the S11 paragraph: in each of the three calculations the two strategies of the
# rule (gee_ar1_fg, mmm_cc) have the fewest expected flags, and the ordering of all strategies is the same.
# Outputs: tables/V5_S8b_sub_s37.csv, tables/V5_S8b_sub_s37_keys.csv, tables/V5_S8b_sub_s37_checks.csv.
# Run from ~/MMM: Rscript --vanilla R/s8b/sub/s8b_table_s37.R
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
SRC <- "tables/V5_S8b_rule_stability.csv"; PROV <- "R/s8b/provenance/lit_pow_rule_stability_v2.R"
CHK <- list(); chk <- function(name, ok, detail = "", fatal = TRUE) { CHK[[length(CHK) + 1]] <<- data.table(check = name, pass = isTRUE(ok), fatal = fatal, detail = detail)
  cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
OUT <- list()
put <- function(key, value, kind, how, mcse = NA_real_, source = SRC) {
  stopifnot(length(value) == 1, !is.na(value))
  OUT[[key]] <<- data.table(key_ = key, value = as.character(value), mcse = if (is.na(mcse)) "" else as.character(mcse), kind = kind, source = source, how = how) }
S <- fread(SRC); T4 <- fread("tables/V5_S7f_table4.csv")
## Table 4 of the main text prints these strategies, in this order (manuscript/v5/MMM_v5.qmd, read-only)
QMD <- readLines("manuscript/v5/MMM_v5.qmd", warn = FALSE)
TL <- grep("^\\| .*f:t4\\|[a-z0-9_]+\\|hi_main", QMD, value = TRUE)
T4S <- unique(sub("^.*?f:t4\\|([a-z0-9_]+)\\|hi_main.*$", "\\1", TL, perl = TRUE))
ROWS <- T4S[T4S %in% S$strategy]
chk("strategies of printed Table 4 found in the stability file", length(ROWS) >= 2 && all(c("gee_ar1_fg", "mmm_cc") %in% ROWS),
    paste0(paste(ROWS, collapse = ","), " | not in file: ", paste(setdiff(T4S, S$strategy), collapse = ",")))
B <- as.integer(sub(".*B <- ([0-9]+)L.*", "\\1", grep("^B <- [0-9]+L", readLines(PROV), value = TRUE)[1]))
chk("bootstrap size read from the provenance script", is.finite(B) && B > 0, as.character(B))
M <- merge(S, T4[, .(strategy, hi_all)], by = "strategy", all.x = TRUE)
chk("observed flags equal Table 4 (hi_obs = hi_all)", M[strategy %in% ROWS, all(hi_obs == hi_all)])
## ---- ranking statement -----------------------------------------------------------------------------------------------
EXP <- c(boot = "mean_flags_boot", plugin = "plugin_exp_flags", eb = "eb_exp_flags")
rk <- function(D, col) D[order(get(col)), strategy]
for (set in c("table4", "file")) { D <- if (set == "table4") S[strategy %in% ROWS] else S
  ords <- lapply(EXP, function(cc) rk(D, cc))
  top2 <- sapply(ords, function(o) setequal(o[1:2], c("gee_ar1_fg", "mmm_cc")))
  chk(sprintf("rule strategies have the fewest expected flags in every calculation (%s rows)", set), all(top2),
      paste(names(EXP), sapply(ords, function(o) paste(o[1:3], collapse = ">")), collapse = " | "), fatal = FALSE)
  same_order <- all(sapply(ords, identical, ords[[1]]))
  put(sprintf("s37|rank_same|%s", set), as.integer(same_order), "int", sprintf("1 = identical ordering by expected flags in the three calculations (%s strategies)", set))
  put(sprintf("s37|rule_top2|%s", set), as.integer(all(top2)), "int", sprintf("1 = gee_ar1_fg and mmm_cc have the fewest expected flags in all three (%s)", set))
  cat(set, "same full ordering in the three calculations:", same_order, "\n") }
pp <- c(boot = "p_pass_boot", plugin = "plugin_p_pass", eb = "eb_p_pass")
po <- lapply(pp, function(cc) S[strategy %in% ROWS][order(-get(cc)), strategy])
put("s37|rank_same_ppass", as.integer(all(sapply(po, identical, po[[1]]))), "int", "1 = identical ordering by P(zero flags) in the three calculations (Table 4 rows)")
## ---- table and keys --------------------------------------------------------------------------------------------------
TT <- S[match(ROWS, strategy), .(strategy, hi_obs, p_pass_boot = 100 * p_pass_boot, mean_flags_boot, q90_flags_boot,
                                 plugin_p_pass = 100 * plugin_p_pass, plugin_exp_flags, eb_p_pass = 100 * eb_p_pass, eb_exp_flags)]
fwrite(TT, "tables/V5_S8b_sub_s37.csv")
for (i in seq_len(nrow(TT))) { r <- TT[i]; s <- r$strategy
  put(sprintf("s37|%s|hi_obs", s), r$hi_obs, "int", "observed design points clearly above 7.5% (Table 4)")
  put(sprintf("s37|%s|p_pass_boot", s), r$p_pass_boot, "sig2", sprintf("%% of %d block-bootstrap resamples with zero flags", B),
      mcse = 100 * sqrt(r$p_pass_boot / 100 * (1 - r$p_pass_boot / 100) / B))
  put(sprintf("s37|%s|mean_flags_boot", s), r$mean_flags_boot, "sig2", "mean number of flags over the resamples")
  put(sprintf("s37|%s|q90_flags_boot", s), r$q90_flags_boot, "int", "90th percentile of the number of flags over the resamples")
  put(sprintf("s37|%s|plugin_p_pass", s), r$plugin_p_pass, "sig2", "binomial plug-in, observed rates: P(zero flags) (%)")
  put(sprintf("s37|%s|plugin_exp_flags", s), r$plugin_exp_flags, "sig2", "binomial plug-in, observed rates: expected flags")
  put(sprintf("s37|%s|eb_p_pass", s), r$eb_p_pass, "sig2", "binomial plug-in, empirical-Bayes rates: P(zero flags) (%)")
  put(sprintf("s37|%s|eb_exp_flags", s), r$eb_exp_flags, "sig2", "binomial plug-in, empirical-Bayes rates: expected flags") }
disp <- function(v) if (v < 0.1) "< 0.1" else if (v >= 10) sprintf("%.0f", v) else sprintf("%.1f", v)
for (i in seq_len(nrow(TT))) for (fl in c("p_pass_boot", "plugin_p_pass", "eb_p_pass"))
  put(sprintf("s37|%s|%s|disp", TT$strategy[i], fl), disp(TT[[fl]][i]), "text", "display of the percentage: '< 0.1' below 0.1, else 1 decimal below 10 and 0 decimals from 10")
put("s37|B", B, "int", "block-bootstrap resamples", source = PROV)
put("s37|nrows", nrow(TT), "int", "strategies of Table 4 in the stability file")
put("s37|fg_pfail_boot", 100 - TT[strategy == "gee_ar1_fg", p_pass_boot], "sig2", "% of resamples with at least one flag, GEE AR(1) FG")
put("s37|mmm_pfail_boot", 100 - TT[strategy == "mmm_cc", p_pass_boot], "sig2", "% of resamples with at least one flag, MMM with check")
KK <- rbindlist(OUT); setnames(KK, "key_", "key"); fwrite(KK, "tables/V5_S8b_sub_s37_keys.csv")
CK <- rbindlist(CHK); fwrite(CK, "tables/V5_S8b_sub_s37_checks.csv")
print(TT, digits = 3); cat("keys", nrow(KK), "| checks failed", sum(!CK$pass), "\n")
stopifnot(all(CK[fatal == TRUE, pass]))
