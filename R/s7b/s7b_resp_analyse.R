# =============================================================================
# s7b_resp_analyse.R -- S7b part 3: measures on the respiratory subsamples with the definitions of rev2/T5_analyse.R and
# rev2/T7_bias.R (usable = finite estimate and SE, |estimate| <= 10, SE < 10; SE ratio = mean SE / SD of estimates;
# bias = mean estimate minus the full-data estimate of the method, the corrections against the GEE).
# Outputs: tables/V4_S7b_resp_usable.csv, tables/V4_S7b_resp_seratio.csv, tables/V4_S7b_resp_bias.csv, tables/V4_S7b_resp_compare.csv
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R"))
A <- fread(file.path(TABDIR, "V4_S7b_resp_fits.csv")); REF <- fread(file.path(TABDIR, "V4_S7b_resp_ref.csv")); DEN <- 500
A[, sep := (is.finite(estimate) & abs(estimate) > 10) | (is.finite(se) & se > 10)]
A[, has_se := is.finite(se) & se > 0]
A[, returned := is.finite(estimate)]
A[, usable := returned & has_se & !sep & abs(estimate) <= 10 & se < 10]
tab <- A[, .(returned = sum(returned), no_se = sum(returned & !has_se), sep_with_se = sum(returned & has_se & sep), usable = sum(usable),
             n_timeout = sum(status == "timeout"), n_error = sum(status == "error"), n_no_term = sum(status == "no_term")), by = .(spec, method, N)]
tab[, `:=`(pct_not_returned = 100 * (DEN - returned) / DEN, pct_no_se_500 = 100 * no_se / DEN, pct_sep_500 = 100 * sep_with_se / DEN,
           pct_usable_500 = 100 * usable / DEN, pct_usable_returned = 100 * usable / returned, mcse_usable_500 = 100 * sqrt((usable / DEN) * (1 - usable / DEN) / DEN))]
setorder(tab, spec, N, -pct_usable_500); fwrite(tab, file.path(TABDIR, "V4_S7b_resp_usable.csv"))
cat("=== usable on 500, by spec / method / N ===\n"); print(dcast(tab, spec + method ~ N, value.var = "pct_usable_500"), digits = 4)
cat("=== usable on returned ===\n"); print(dcast(tab, spec + method ~ N, value.var = "pct_usable_returned"), digits = 4)
U <- A[usable == TRUE]
R <- U[, .(n = .N, mean_se = mean(se), sd_est = sd(estimate)), by = .(spec, method, N)][, ratio := mean_se / sd_est]
R[, mcse_ratio := ratio * sqrt(1 / (2 * (n - 1)))]      # from the SD's MCSE only (Morris Tab. 6, empSE), the mean SE's own error is far smaller
fwrite(R, file.path(TABDIR, "V4_S7b_resp_seratio.csv"))
cat("=== mean SE / SD of estimates ===\n"); print(dcast(R, spec + method ~ N, value.var = "ratio"), digits = 3)
RR <- REF[, .(spec, method, ref = estimate)]; rg <- REF[method == "GEE", .(spec, ref_gee = estimate)]
B <- U[, .(n = .N, mean_est = mean(estimate), sd_est = sd(estimate)), by = .(spec, method, N)]
B <- merge(B, RR, by = c("spec", "method"), all.x = TRUE); B <- merge(B, rg, by = "spec", all.x = TRUE)
B[method %in% c("GEE-FG", "GEE-KC", "GEE-MD"), ref := ref_gee]
B[, `:=`(bias = mean_est - ref, mcse_bias = sd_est / sqrt(n), bias_over_sd = (mean_est - ref) / sd_est)]
setorder(B, spec, method, N); fwrite(B, file.path(TABDIR, "V4_S7b_resp_bias.csv"))
cat("=== bias against the full-data estimate ===\n"); print(dcast(B, spec + method ~ N, value.var = "bias"), digits = 3)
# ---- side by side with the neonatal cohort (same definitions) -----------------------------------------------------------
L5 <- fread(file.path(REPO, "rev2/T5_tab316.csv")); L5r <- fread(file.path(REPO, "rev2/T5_seratio.csv")); L7 <- fread(file.path(REPO, "rev2/T7_bias.csv"))
CMP <- rbind(merge(tab[spec %in% c("pieno", "ridotto"), .(spec, method, N, resp_usable_500 = pct_usable_500, resp_usable_returned = pct_usable_returned)],
                   L5[spec %in% c("pieno", "ridotto"), .(spec, method, N, lus_usable_500 = pct_usable_500, lus_usable_returned = pct_usable_returned)], by = c("spec", "method", "N"), all = TRUE))
CMP <- merge(CMP, R[, .(spec, method, N, resp_seratio = ratio)], by = c("spec", "method", "N"), all.x = TRUE)
CMP <- merge(CMP, L5r[, .(spec, method, N, lus_seratio = ratio)], by = c("spec", "method", "N"), all.x = TRUE)
CMP <- merge(CMP, B[, .(spec, method, N, resp_bias = bias, resp_bias_over_sd = bias_over_sd)], by = c("spec", "method", "N"), all.x = TRUE)
CMP <- merge(CMP, L7[, .(spec, method, N, lus_bias = bias)], by = c("spec", "method", "N"), all.x = TRUE)
setorder(CMP, spec, method, N); fwrite(CMP, file.path(TABDIR, "V4_S7b_resp_compare.csv"))
print(CMP[spec == "pieno" & N %in% c(10, 50)], digits = 3)
