# =============================================================================
# s7b_preds.R -- S7b: judgement of the pre-registered predictions PC1-PC12 (docs/steps/S7b.md) from the step tables,
# and the finite-population view of the SE ratio in both applications (an anomaly met during the work, see the report).
# Outputs: tables/V4_S7b_preds.csv, tables/V4_S7b_seratio_fpc.csv
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R"))
judge <- function(est, mc, lo, hi) { a <- est - 2 * mc; b <- est + 2 * mc
  fifelse(a >= lo & b <= hi, "confermata", fifelse(b < lo | a > hi, "smentita", "indeterminata")) }
SM <- fread(file.path(TABDIR, "V4_S7b_rule_summary.csv")); SH <- fread(file.path(TABDIR, "V4_S7b_df_shift.csv"))
NS <- fread(file.path(TABDIR, "V4_S7b_null_summary.csv")); NSd <- fread(file.path(TABDIR, "V4_S7b_null_shift.csv"))
DG <- fread(file.path(TABDIR, "V4_S7b_resp_degeneracy.csv")); US <- fread(file.path(TABDIR, "V4_S7b_resp_usable.csv"))
SR <- fread(file.path(TABDIR, "V4_S7b_resp_seratio.csv")); BI <- fread(file.path(TABDIR, "V4_S7b_resp_bias.csv"))
P <- list()
add <- function(code, what, est, mc, lo, hi, source, note = "") P[[length(P) + 1L]] <<- data.table(code = code, quantity = what, observed = est, mcse = mc, rel_lo = lo, rel_hi = hi, verdict = judge(est, mc, lo, hi), source = source, note = note)
# PC1, PC2
x <- SM[src == "tutti" & strategy == "mmm|t2"]
add("PC1", "design points with MMM (check) + t at N - 2 clearly above 7.5%", x$n_hi, 0, 0, 0, "V4_S7b_rule_summary.csv, strategy mmm|t2, src tutti, n_hi")
rec <- SM[src == "tutti" & ref == "t2" & recommended == TRUE, strategy]
add("PC2", "strategy chosen by the rule with N - 2 (1 = mmm|t2, 0 = other)", as.numeric(identical(rec, "mmm|t2")), 0, 1, 1, "V4_S7b_rule_summary.csv, ref t2, recommended", paste("chosen:", paste(rec, collapse = ",")))
# PC3: every Wald strategy, all design points
for (m in SH[group == "all design points", method]) { s <- SH[group == "all design points" & method == m]; add("PC3", sprintf("mean shift of the type I error rate, N - 2 minus N - p, all design points, %s (points)", m), s$shift, s$mcse, 0, 1, "V4_S7b_df_shift.csv, group all design points") }
for (m in c("mmm", "gee_fg")) { s <- SH[group == "N = 10" & method == m]; add("PC4", sprintf("mean shift at N = 10, %s (points)", m), s$shift, s$mcse, 0, 2.5, "V4_S7b_df_shift.csv, group N = 10") }
x <- SM[src == "tutti" & strategy == "gee_fg|t2"]
add("PC5", "design points with Fay-Graubard + t at N - 2 clearly above 7.5%", x$n_hi, 0, 2, 10, "V4_S7b_rule_summary.csv, strategy gee_fg|t2, src tutti, n_hi")
# PC6, PC7: max excess (deterministic given the cell SE; the MCSE of the empirical SE is propagated: excess scales with |delta|/s, MCSE of s = s/sqrt(2(n-1)))
X <- NSd[ref == "n"]; X[, mc_exc := abs(excess) * (1 / sqrt(2 * (n - 1)))]
i <- X[which.max(excess)]; add("PC6", "largest expected excess rejection over all GLMM-native design points and marginal methods, normal reference (points)", i$excess, i$mc_exc, 0, 1.5, "V4_S7b_null_shift.csv, ref n, max excess", sprintf("%s N=%d T=%d %s; empirical SE %.4f, rate %.1f", i$scen, i$N, i$ntime, i$method, i$empse, i$rate))
j <- X[method != "shrink"][which.max(excess)]; add("PC6b", "the same, shrinkage excluded (not pre-registered, explanation of PC6)", j$excess, j$mc_exc, 0, 1.5, "V4_S7b_null_shift.csv, ref n, max excess, method != shrink", sprintf("%s N=%d T=%d %s; empirical SE %.4f, rate %.1f", j$scen, j$N, j$ntime, j$method, j$empse, j$rate))
k <- X[group != "M6-M7"][which.max(excess)]; add("PC7", "largest expected excess rejection, main scenarios M1-M5 and M1c, normal reference (points)", k$excess, k$mc_exc, 0, 0.5, "V4_S7b_null_shift.csv, ref n, group M1-M5 and M1c", sprintf("%s N=%d T=%d %s; empirical SE %.4f", k$scen, k$N, k$ntime, k$method, k$empse))
for (m in c("mmm", "gee_fg")) { s <- NS[group == "all" & ref == "t" & method == m]; add("PC8", sprintf("hi flags on GLMM-native that change after subtracting the excess, %s + t", m), s$n_hi_flip, 0, 0, 0, "V4_S7b_null_summary.csv, group all, ref t") }
# PC9-PC12
p9 <- 100 * DG[N == 10, mean(any_one)]; add("PC9", "share of draws of 10 patients with a binary covariate at one level (%)", p9, 100 * sqrt(p9 / 100 * (1 - p9 / 100) / 500), 0, 25, "V4_S7b_resp_degeneracy.csv, N = 10, any_one")
u <- US[spec == "pieno" & method == "MMM" & N == 10]; p10 <- u$pct_usable_returned; add("PC10", "usable on returned, 10 patients, full specification, MMM (%)", p10, 100 * sqrt(p10 / 100 * (1 - p10 / 100) / u$returned), 75, 100, "V4_S7b_resp_usable.csv, pieno MMM N=10, pct_usable_returned", sprintf("returned %d, usable %d of 500", u$returned, u$usable))
s <- SR[spec == "pieno" & method == "MMM" & N == 50]; add("PC11", "mean model SE / SD of estimates, MMM, full specification, 50 patients", s$ratio, s$mcse_ratio, 0.91, 1.14, "V4_S7b_resp_seratio.csv, pieno MMM N=50, ratio")
b <- BI[spec == "pieno" & method == "MMM" & N == 50]; add("PC12", "bias / SD of estimates, MMM, full specification, 50 patients", b$bias_over_sd, 1 / sqrt(b$n), -0.46, 0.46, "V4_S7b_resp_bias.csv, pieno MMM N=50, bias_over_sd")
PR <- rbindlist(P); fwrite(PR, file.path(TABDIR, "V4_S7b_preds.csv")); print(PR[, .(code, observed = signif(observed, 3), mcse = signif(mcse, 2), rel_lo, rel_hi, verdict)])
# ---- finite-population view of the SE ratio: subsamples without replacement from P patients have SD shrunk by sqrt(1 - N/P) ----
L5r <- fread(file.path(REPO, "rev2/T5_seratio.csv"))[spec %in% c("pieno", "ridotto")][, `:=`(application = "neonatal cohort", P = 207L)]
R2 <- SR[, .(spec, method, N, n, mean_se, sd_est, ratio)][, `:=`(application = "respiratory", P = 111L)]
F <- rbind(L5r[, .(application, spec, method, N, P, n, mean_se, sd_est, ratio)], R2[, .(application, spec, method, N, P, n, mean_se, sd_est, ratio)])
F[, `:=`(fpc = sqrt(1 - N / P), ratio_fpc = ratio * sqrt(1 - N / P), mcse_ratio = ratio / sqrt(2 * (n - 1)))]
fwrite(F, file.path(TABDIR, "V4_S7b_seratio_fpc.csv"))
cat("\nSE ratio with the finite-population factor, full specification, MMM and GLMM:\n")
print(dcast(F[spec == "pieno" & method %in% c("MMM", "GLMM", "GEE-FG")], application + method ~ N, value.var = "ratio_fpc"), digits = 3)
print(dcast(F[spec == "pieno" & method %in% c("MMM", "GLMM", "GEE-FG")], application + method ~ N, value.var = "ratio"), digits = 3)
