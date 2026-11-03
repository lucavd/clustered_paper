#!/usr/bin/env Rscript
# =============================================================================
# s7f_kc_clinical_tables.R -- S7f: tables of the two clinical applications with the Kauermann-Carroll rows from fewclust.
# The rows of every other method are the stored ones; the measures use the definitions of rev2/T5_analyse.R,
# rev2/T7_bias.R, R/s7b/s7b_resp_analyse.R and R/s7b/s7b_preds.R (finite-population factor), unchanged.
# Checks (written to tables/V5_S7f_kc_clinical_checks.csv):
#   (a) the GEE estimate refitted here equals the stored one in every draw where both exist;
#   (b) fewclust Mancl-DeRouen and Fay-Graubard against the stored geesmv values: share of draws within 1%;
#   (c) every non-KC row of the rebuilt tables equals the stored table.
# Inputs read-only: rev2/T3_subsample_rerun.csv, rev2/T3_ref.csv, rev2/T5_tab316.csv, rev2/T5_seratio.csv, rev2/T7_bias.csv,
#   tables/V4_S7b_resp_fits.csv, tables/V4_S7b_resp_ref.csv, tables/V4_S7b_resp_usable.csv, V4_S7b_resp_bias.csv,
#   tables/V4_S7b_seratio_fpc.csv, tables/V5_S7f_lus_fc.csv, tables/V5_S7f_resp_fc.csv.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
AGREE <- list()
CHK <- list(); chk <- function(app, what, value, pass) CHK[[length(CHK) + 1L]] <<- data.table(app = app, check = what, value = value, pass = pass)
same_tab <- function(a, b, keys) {   # numeric columns equal (or both NA) on the common rows
  m <- merge(a, b, by = keys, suffixes = c(".a", ".b")); cols <- intersect(setdiff(names(a), keys), setdiff(names(b), keys))
  bad <- 0L
  for (cc in cols) { x <- m[[paste0(cc, ".a")]]; y <- m[[paste0(cc, ".b")]]
    if (is.numeric(x) && is.numeric(y)) bad <- bad + sum(!((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) < 1e-9)))
    else bad <- bad + sum(!((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & as.character(x) == as.character(y)))) }
  c(rows = nrow(m), bad = bad)
}
swap_kc <- function(A, REF, FC, app, extra_keys) {
  K <- c("spec", "N", "rep")
  # (a) GEE estimate
  g <- merge(A[method == "GEE", c(K, "estimate"), with = FALSE], unique(FC[rep > 0 & corr == "kc", c(K, "gee_est"), with = FALSE]), by = K)
  both <- g[is.finite(estimate) & is.finite(gee_est)]
  chk(app, "GEE estimate refitted = stored (draws compared; max abs diff)", sprintf("%d; %.3g", nrow(both), if (nrow(both)) max(abs(both$estimate - both$gee_est)) else NA),
      nrow(both) > 0 && max(abs(both$estimate - both$gee_est)) < 1e-8)
  # (b) fewclust against geesmv for MD and FG
  for (cr in c("md", "fg")) {
    m <- merge(A[method == paste0("GEE-", toupper(cr)) & status == "ok", c(K, "se"), with = FALSE], FC[rep > 0 & corr == cr & status == "ok", c(K, "se"), with = FALSE], by = K, suffixes = c(".smv", ".fc"))
    chk(app, sprintf("%s: fewclust within 1%% of geesmv (draws compared; share %%)", toupper(cr)), sprintf("%d; %.2f", nrow(m), 100 * mean(abs(m$se.fc / m$se.smv - 1) < 0.01)), TRUE)
    AGREE[[length(AGREE) + 1L]] <<- data.table(app = app, corr = cr, n = nrow(m), share_within_1pct = 100 * mean(abs(m$se.fc / m$se.smv - 1) < 0.01))
  }
  kc <- FC[rep > 0 & corr == "kc"]
  B <- A[method == "GEE-KC"]; B[, rid := .I]
  B <- merge(B, kc[, c(K, "se", "status", "gee_est"), with = FALSE], by = K, all.x = TRUE, suffixes = c("", ".fc"))
  chk(app, "GEE-KC rows without a fewclust row", B[is.na(status.fc), .N], B[is.na(status.fc), .N] == 0L)
  B[, `:=`(se = fifelse(status.fc == "ok", se.fc, NA_real_), status = status.fc)]
  B[, c("se.fc", "status.fc", "gee_est", "rid") := NULL]
  A2 <- rbind(A[method != "GEE-KC"], B[, names(A), with = FALSE])
  R2 <- copy(REF); r <- FC[rep == 0 & corr == "kc"]
  for (s in r$spec) R2[method == "GEE-KC" & spec == s, `:=`(se = r[spec == s, if (status == "ok") se else NA_real_], status = r[spec == s, status])]
  list(A = A2, REF = R2, kc = kc)
}
## ---- neonatal cohort -------------------------------------------------------------------------------------------------
A <- fread("rev2/T3_subsample_rerun.csv"); REF <- fread("rev2/T3_ref.csv"); FC <- fread("tables/V5_S7f_lus_fc.csv")
S <- swap_kc(A, REF, FC, "neonatal cohort"); A <- S$A; REF <- S$REF; DEN <- 500
A[, sep := is.finite(estimate) & abs(estimate) > 10 | (is.finite(se) & se > 10)]
A[, has_se := is.finite(se) & se > 0]
A[, returned := is.finite(estimate)]
A[, usable := returned & has_se & !sep & is.finite(estimate) & abs(estimate) <= 10 & se < 10]
tab <- A[, .(returned = sum(returned), no_se = sum(returned & !has_se), sep_with_se = sum(returned & has_se & sep), usable = sum(usable),
             conv_only = sum(returned)), by = .(spec, method, N)]
tab[, `:=`(pct_not_returned = 100 * (DEN - returned) / DEN, pct_no_se_500 = 100 * no_se / DEN, pct_sep_500 = 100 * sep_with_se / DEN,
           pct_usable_500 = 100 * usable / DEN, pct_no_se_returned = 100 * no_se / returned, pct_usable_returned = 100 * usable / returned,
           pct_convonly_500 = 100 * conv_only / DEN)]
tab[, somma := pct_not_returned + pct_no_se_500 + pct_sep_500 + pct_usable_500]
fwrite(tab, "tables/V5_S7f_lus_tab316.csv")
U <- A[usable == TRUE]
Rr <- U[, .(n = .N, mean_se = mean(se), sd_est = sd(estimate)), by = .(spec, method, N)][, ratio := mean_se / sd_est]
fwrite(Rr, "tables/V5_S7f_lus_seratio.csv")
Bl <- U[, .(n = .N, mean_est = mean(estimate)), by = .(spec, method, N)]
Bl <- merge(Bl, REF[, .(spec, method, ref = estimate)], by = c("spec", "method"), all.x = TRUE)
Bl <- merge(Bl, REF[method == "MMM", .(spec, ref_mmm = estimate)], by = "spec", all.x = TRUE)
Bl <- merge(Bl, REF[method == "GEE", .(spec, ref_gee = estimate)], by = "spec", all.x = TRUE)
Bl[method %in% c("GEE-FG", "GEE-KC", "GEE-MD"), ref := ref_gee]
Bl[, `:=`(bias = mean_est - ref, bias_vs_mmm = mean_est - ref_mmm)]
fwrite(Bl, "tables/V5_S7f_lus_bias.csv"); fwrite(REF, "tables/V5_S7f_lus_ref.csv")
x <- same_tab(fread("rev2/T5_tab316.csv")[method != "GEE-KC"], tab[method != "GEE-KC"], c("spec", "method", "N"))
chk("neonatal cohort", "tab316: non-KC rows equal to rev2/T5_tab316.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
x <- same_tab(fread("rev2/T5_seratio.csv")[method != "GEE-KC"], Rr[method != "GEE-KC"], c("spec", "method", "N"))
chk("neonatal cohort", "seratio: non-KC rows equal to rev2/T5_seratio.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
x <- same_tab(fread("rev2/T7_bias.csv")[method != "GEE-KC"], Bl[method != "GEE-KC"], c("spec", "method", "N"))
chk("neonatal cohort", "bias: non-KC rows equal to rev2/T7_bias.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
## ---- respiratory trial -------------------------------------------------------------------------------------------------
A <- fread("tables/V4_S7b_resp_fits.csv"); REF <- fread("tables/V4_S7b_resp_ref.csv"); FC <- fread("tables/V5_S7f_resp_fc.csv")
S <- swap_kc(A, REF, FC, "respiratory"); A <- S$A; REF <- S$REF
A[, sep := (is.finite(estimate) & abs(estimate) > 10) | (is.finite(se) & se > 10)]
A[, has_se := is.finite(se) & se > 0]
A[, returned := is.finite(estimate)]
A[, usable := returned & has_se & !sep & abs(estimate) <= 10 & se < 10]
tb <- A[, .(returned = sum(returned), no_se = sum(returned & !has_se), sep_with_se = sum(returned & has_se & sep), usable = sum(usable),
            n_timeout = sum(status == "timeout"), n_error = sum(status == "error"), n_no_term = sum(status == "no_term")), by = .(spec, method, N)]
tb[, `:=`(pct_not_returned = 100 * (DEN - returned) / DEN, pct_no_se_500 = 100 * no_se / DEN, pct_sep_500 = 100 * sep_with_se / DEN,
          pct_usable_500 = 100 * usable / DEN, pct_usable_returned = 100 * usable / returned, mcse_usable_500 = 100 * sqrt((usable / DEN) * (1 - usable / DEN) / DEN))]
setorder(tb, spec, N, -pct_usable_500); fwrite(tb, "tables/V5_S7f_resp_usable.csv")
U <- A[usable == TRUE]
Rs <- U[, .(n = .N, mean_se = mean(se), sd_est = sd(estimate)), by = .(spec, method, N)][, ratio := mean_se / sd_est]
Rs[, mcse_ratio := ratio * sqrt(1 / (2 * (n - 1)))]
Bs <- U[, .(n = .N, mean_est = mean(estimate), sd_est = sd(estimate)), by = .(spec, method, N)]
Bs <- merge(Bs, REF[, .(spec, method, ref = estimate)], by = c("spec", "method"), all.x = TRUE)
Bs <- merge(Bs, REF[method == "GEE", .(spec, ref_gee = estimate)], by = "spec", all.x = TRUE)
Bs[method %in% c("GEE-FG", "GEE-KC", "GEE-MD"), ref := ref_gee]
Bs[, `:=`(bias = mean_est - ref, mcse_bias = sd_est / sqrt(n), bias_over_sd = (mean_est - ref) / sd_est)]
setorder(Bs, spec, method, N); fwrite(Bs, "tables/V5_S7f_resp_bias.csv"); fwrite(REF, "tables/V5_S7f_resp_ref.csv")
x <- same_tab(fread("tables/V4_S7b_resp_usable.csv")[method != "GEE-KC"], tb[method != "GEE-KC"], c("spec", "method", "N"))
chk("respiratory", "usable: non-KC rows equal to V4_S7b_resp_usable.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
x <- same_tab(fread("tables/V4_S7b_resp_bias.csv")[method != "GEE-KC"], Bs[method != "GEE-KC"], c("spec", "method", "N"))
chk("respiratory", "bias: non-KC rows equal to V4_S7b_resp_bias.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
x <- same_tab(fread("tables/V4_S7b_resp_ref.csv")[method != "GEE-KC"], REF[method != "GEE-KC"], c("spec", "method"))
chk("respiratory", "full-data estimates: non-KC rows equal to V4_S7b_resp_ref.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
## ---- SE ratio with the finite-population factor, both applications (definition of R/s7b/s7b_preds.R) -------------------
L5r <- Rr[spec %in% c("pieno", "ridotto")][, `:=`(application = "neonatal cohort", P = 207L)]
R2 <- Rs[, .(spec, method, N, n, mean_se, sd_est, ratio)][, `:=`(application = "respiratory", P = 111L)]
F <- rbind(L5r[, .(application, spec, method, N, P, n, mean_se, sd_est, ratio)], R2[, .(application, spec, method, N, P, n, mean_se, sd_est, ratio)])
F[, `:=`(fpc = sqrt(1 - N / P), ratio_fpc = ratio * sqrt(1 - N / P), mcse_ratio = ratio / sqrt(2 * (n - 1)))]
fwrite(F, "tables/V5_S7f_seratio_fpc.csv")
x <- same_tab(fread("tables/V4_S7b_seratio_fpc.csv")[method != "GEE-KC"], F[method != "GEE-KC"], c("application", "spec", "method", "N"))
chk("both", "SE ratio with fpc: non-KC rows equal to V4_S7b_seratio_fpc.csv (rows; discordant cells)", sprintf("%d; %d", x[1], x[2]), x[2] == 0)
CK <- rbindlist(CHK); fwrite(CK, "tables/V5_S7f_kc_clinical_checks.csv"); print(CK)
fwrite(rbindlist(AGREE), "tables/V5_S7f_kc_clinical_agree.csv")
cat("\nKauermann-Carroll, full model, by N (fewclust):\n")
print(tab[spec == "pieno" & method %in% c("GEE-KC", "GEE-FG", "GEE-MD"), .(method, N, pct_usable_500, pct_no_se_500)][order(method, N)])
print(tb[spec == "pieno" & method %in% c("GEE-KC", "GEE-FG", "GEE-MD"), .(method, N, pct_usable_500, pct_no_se_500)][order(method, N)])
print(F[spec == "pieno" & method %in% c("GEE-KC", "GEE-FG"), .(application, method, N, ratio_fpc)][order(application, method, N)])
print(REF[method %in% c("GEE", "GEE-KC", "GEE-FG")])
stopifnot(all(CK$pass))
