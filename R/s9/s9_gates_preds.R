# =============================================================================
# s9_gates_preds.R -- S9: gates GS2, GS3 and predictions PS1-PS8 (pre-registered in docs/steps/S9.md), and the diff of the
# number registry before and after the step. "Before" = results_v4/S9/pre/tables (extracted from tables_pre_S9.tgz).
# Outputs: tables/V6_S9_GS2.csv, V6_S9_GS3.csv, V6_S9_predictions.csv, V6_S9_registry_diff.csv. Read-only on inputs.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
PRE <- "results_v4/S9/pre/tables"; NEW <- "tables"
rd <- function(d, f) fread(file.path(d, f))
## ---- GS2 -------------------------------------------------------------------------------------------------------------
same_tab <- function(A, B, keys, keep) {
  A <- A[keep(A)]; B <- B[keep(B)]
  cols <- setdiff(intersect(names(A), names(B)), keys)
  M <- merge(A, B, by = keys, suffixes = c(".new", ".old"))
  bad <- 0L
  for (cc in cols) { a <- M[[paste0(cc, ".new")]]; b <- M[[paste0(cc, ".old")]]
    if (is.numeric(a) && is.numeric(b)) ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12)
    else ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
    bad <- bad + sum(!ok) }
  data.table(rows_new = nrow(A), rows_old = nrow(B), rows_matched = nrow(M), cells_differing = bad)
}
U7d <- c("gee_exch_cr2", "gee_ind_cr2", "flip_strict", "study_mmm", "study_gee_md", "study_gee_kc", "study_gee_fg")
U7e <- c("mmm_cc", "mmm_study", "glmer", "shrink", "glmer_lrt", "gee_exch_cr2", "gee_ind_cr2", "flip_strict")
URD <- c("mmm_cc", "glmer", "shrink", "glmer_lrt", "flip", "gee_exch_cr2", "gee_ind_cr2")
G <- list(
  cbind(table = "V4_S7d_dp_rates (S7d)", same_tab(rd(NEW, "V4_S7d_dp_rates.csv"), rd(PRE, "V4_S7d_dp_rates.csv"), c("dgp", "scen", "N", "ntime", "strategy", "ref"), function(X) X$strategy %in% U7d)),
  cbind(table = "V4_S7e_dp_rates (S7e)", same_tab(rd(NEW, "V4_S7e_dp_rates.csv"), rd(PRE, "V4_S7e_dp_rates.csv"), c("dgp", "scen", "N", "ntime", "strategy", "ref"), function(X) X$strategy %in% U7e)),
  cbind(table = "V4_S7b_rule_dp (S7b)", same_tab(rd(NEW, "V4_S7b_rule_dp.csv"), rd(PRE, "V4_S7b_rule_dp.csv"), c("dgp", "scen", "N", "ntime", "strategy"), function(X) X$method != "gee")),
  cbind(table = "V4_S6b_ns_rule_safety_dp (M6-M7)", same_tab(rd(NEW, "V4_S6b_ns_rule_safety_dp.csv"), rd(PRE, "V4_S6b_ns_rule_safety_dp.csv"), c("dgp", "scen", "N", "ntime", "strategy"), function(X) !startsWith(X$strategy, "gee|"))),
  cbind(table = "variant_nc vs kcfc V4_S6_rule_safety_dp", same_tab(rd(NEW, "variant_nc/V4_S6_rule_safety_dp.csv"), rd(PRE, "variant_kcfc/V4_S6_rule_safety_dp.csv"), c("dgp", "scen", "N", "ntime", "strategy"), function(X) !startsWith(X$strategy, "gee|") & !startsWith(X$strategy, "gee_kc|"))),
  cbind(table = "variant_nc vs kcfc V4_perf_all (own set)", same_tab(rd(NEW, "variant_nc/V4_perf_all.csv.gz"), rd(PRE, "variant_kcfc/V4_perf_all.csv.gz"), c("dgp", "scen", "model", "N", "ntime", "label", "arm", "method", "usab", "set"), function(X) !X$method %in% c("gee", "gee_kc") & X$set == "own")),
  cbind(table = "V5_S7f_rule_dp (Table 5 per point)", same_tab(rd(NEW, "V5_S7f_rule_dp.csv"), rd(PRE, "V5_S7f_rule_dp.csv"), c("dgp", "scen", "N", "ntime", "strategy"), function(X) X$strategy %in% URD)),
  cbind(table = "V4_S7d2_dp (sign-flip test columns)", same_tab(rd(NEW, "V4_S7d2_dp.csv")[, .(dgp, scen, N, ntime, label, n_ut, cov_t, pu_flip_s7d, pu_flip_guard, rate_flip_guard, mcse_cov_t, lo_t, hi_t)],
                                                          rd(PRE, "V4_S7d2_dp.csv")[, .(dgp, scen, N, ntime, label, n_ut, cov_t, pu_flip_s7d, pu_flip_guard, rate_flip_guard, mcse_cov_t, lo_t, hi_t)], c("dgp", "scen", "N", "ntime", "label"), function(X) rep(TRUE, nrow(X)))))
GS2 <- rbindlist(G); GS2[, pass := rows_new == rows_old & rows_matched == rows_new & cells_differing == 0L]
print(GS2); fwrite(GS2, file.path(NEW, "V6_S9_GS2.csv"))
## ---- GS3 -------------------------------------------------------------------------------------------------------------
CK <- rd(NEW, "V6_S9_const_check_table5.csv"); T5 <- rd(NEW, "V5_S7f_table4.csv")
M <- merge(CK[, .(strategy, hi_corr, lo_corr, pu_corr, pu_small_corr, upow_corr)], T5[, .(strategy, hi_all, lo_all, pu, pu_small, upow)], by = "strategy")
M[, pass := hi_corr == hi_all & lo_corr == lo_all & abs(pu_corr - pu) < 1e-9 & abs(pu_small_corr - pu_small) < 1e-9 & abs(upow_corr - upow) < 1e-9]
print(M, digits = 6); fwrite(M, file.path(NEW, "V6_S9_GS3.csv")); cat("GS3:", M[, sum(pass)], "of", nrow(M), "\n")
## ---- registry diff ----------------------------------------------------------------------------------------------------
RO <- rd(PRE, "V5_S7f_numbers.csv"); RN <- rd(NEW, "V5_S7f_numbers.csv")
num <- function(x) suppressWarnings(as.numeric(x))
RD <- merge(RO[, .(id, value_old = num(value), text_old = text, scale_old = num(scale))], RN[, .(id, file, value_new = num(value), text_new = text, scale = num(scale))], by = "id", all = TRUE)
RD[, status := fifelse(is.na(text_old), "new id", fifelse(is.na(text_new), "id removed", fifelse(text_old == text_new, fifelse(abs(value_new - value_old) < 1e-12 | (is.na(value_new) & is.na(value_old)), "same", "same text, value moved"), "text changed")))]
fwrite(RD, file.path(NEW, "V6_S9_registry_diff.csv")); print(RD[, .N, by = status])
## ---- predictions -----------------------------------------------------------------------------------------------------
dv <- function(ids) { X <- RD[id %in% ids]; stopifnot(nrow(X) == length(ids)); X[, .(id, old = value_old * scale, new = value_new * scale, d = (value_new - value_old) * scale)] }
verdict <- function(lo, hi, a, b) fifelse(lo >= a & hi <= b, "confermata", fifelse(hi < a | lo > b, "smentita", "indeterminata"))
P <- list()
EG <- expand.grid(g = c("GB", "NAT", "AR1"), m = c("gee", "gee_kc"), w = c("rate", "mad"), stringsAsFactors = FALSE)
x <- dv(paste0("v:", EG$w, "|", EG$g, "|", EG$m, "@tab"))
P[["PS1"]] <- x[, .(code = "PS1", what = "Table 3 rate and MAD, GEE uncorrected and KC (12 values)", d_min = min(d), d_max = max(d), a = -1, b = 1)]
x <- dv(as.vector(outer(c("GB", "NAT", "AR1"), c("gee", "gee_kc"), function(g, m) paste0("v:wc|", g, "|", m, "@tab"))))
P[["PS2"]] <- x[, .(code = "PS2", what = "Table 3 worst case, GEE uncorrected and KC (6 values)", d_min = min(d), d_max = max(d), a = -2.5, b = 2.5)]
x <- dv(as.vector(outer(c("RI", "RS"), c(2, 3, 4, 6, 10), function(r, t) paste0("v:puT|", r, "|gee|", t, "@tab"))))
P[["PS3"]] <- x[, .(code = "PS3", what = "Table 4 usable share, GEE uncorrected (10 cells)", d_min = min(d), d_max = max(d), a = -2.5, b = 2.5)]
P[["PS3b"]] <- x[grepl("\\|(4|6|10)@tab$", id), .(code = "PS3 (T >= 4)", what = "same, T from 4 (expected -0.3 to 0)", d_min = min(d), d_max = max(d), a = -2.5, b = 2.5)]
TCo <- rd(PRE, "variant_kcfc/V4_S6_table_cells.csv"); TCn <- rd(NEW, "variant_nc/V4_S6_table_cells.csv")
tc <- merge(TCn[grepl("^T10", table_id)], TCo[grepl("^T10", table_id)], by = intersect(c("table_id", "variant", "i", "j"), names(TCn)), suffixes = c(".n", ".o"))
P[["PS4"]] <- tc[, .(code = "PS4", what = sprintf("Table S10 power cells (%d)", .N), d_min = min(est.n - est.o, na.rm = TRUE), d_max = max(est.n - est.o, na.rm = TRUE), a = -2, b = 2)]
Ho <- rd(PRE, "V5_S7f_h2h_boot.csv"); Hn <- rd(NEW, "V5_S7f_h2h_boot.csv")
H <- merge(Hn, Ho, by = c("region", "q"), suffixes = c(".n", ".o"))
Hd <- H[grepl("mad_gee_minus_mmm", q)]
P[["PS5"]] <- Hd[, .(code = "PS5", what = sprintf("direct comparison, MAD difference GEE - MMM, %d regions", .N), d_min = min(est.n - est.o), d_max = max(est.n - est.o), a = -1, b = 1)]
cls <- H[!is.na(class.n) | !is.na(class.o), .(region, q, class.o, class.n)][class.o != class.n | xor(is.na(class.o), is.na(class.n))]
fwrite(H, file.path(NEW, "V6_S9_h2h_before_after.csv"))
Bo <- rd(PRE, "V4_S7e_bias_summary.csv"); Bn <- rd(NEW, "V4_S7e_bias_summary.csv"); B <- merge(Bn, Bo, by = "strategy", suffixes = c(".n", ".o"))[startsWith(strategy, "gee")]
P[["PS6"]] <- B[, .(code = "PS6", what = sprintf("copula, mean coverage of the GEE strategies (t and normal, %d strategies)", .N), d_min = min(c(mean_cov_t.n - mean_cov_t.o, mean_cov_n.n - mean_cov_n.o)), d_max = max(c(mean_cov_t.n - mean_cov_t.o, mean_cov_n.n - mean_cov_n.o)), a = -1, b = 1)]
R2 <- rd(NEW, "V5_S7f_rule_t2.csv")[df == "N - 2"]
P[["PS7"]] <- data.table(code = "PS7", what = sprintf("N - 2: GEE AR(1) FG clearly above = %d; strategy chosen = %s", R2[strategy == "gee_ar1_fg", n_hi], paste(R2[chosen == TRUE, strategy], collapse = " ")),
                         d_min = as.numeric(R2[strategy == "gee_ar1_fg", n_hi] >= 2 & identical(R2[chosen == TRUE, strategy], "mmm_cc")), d_max = NA_real_, a = 1, b = 1)
RS <- rd(NEW, "V5_S8b_rule_stability.csv")
w <- c(boot = RS[which.min(mean_flags_boot), strategy], plugin = RS[which.min(plugin_exp_flags), strategy], eb = RS[which.min(eb_exp_flags), strategy])
P[["PS8"]] <- data.table(code = "PS8", what = paste("fewest expected flags:", paste(names(w), w, sep = "=", collapse = ", ")), d_min = as.numeric(all(w == "mmm_cc")), d_max = NA_real_, a = 1, b = 1)
PR <- rbindlist(P, fill = TRUE)
PR[, verdict := fifelse(code %in% c("PS7", "PS8"), fifelse(d_min == 1, "confermata", "smentita"), verdict(d_min, d_max, a, b))]
print(PR, digits = 3); fwrite(PR, file.path(NEW, "V6_S9_predictions.csv"))
cat("Figure 1a classes that changed:\n"); print(cls)
fwrite(RS, file.path(NEW, "V6_S9_rule_stability_view.csv"))
