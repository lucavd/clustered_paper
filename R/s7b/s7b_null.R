# =============================================================================
# s7b_null.R -- S7b part 2: expected excess rejection due to the non-zero marginal null under GLMM-native.
# For each GLMM-native design point, H0 arm and marginal method: delta = marginal truth under H0 of the cell
# (tables/V4_S4_truth_by_cell.csv, tables/V4_S6b_ns_truth.csv), s = empirical SE of the usable estimates,
# excess = P(|Z + delta/s| > c) - P(|Z| > c) with c the normal (ref n) or t_{N-p} (ref t) critical value.
# The observed rate minus the excess is re-flagged with the Bradley limits. Gate GC4: Monte Carlo counter-check.
# Input : tables/V4_S7b_dp_rates.csv, tables/V4_S7b_rule_dp.csv (s7b_df.R)
# Outputs: tables/V4_S7b_null_shift.csv, tables/V4_S7b_null_summary.csv, tables/V4_S7b_null_flips.csv, GC4 row in tables/V4_S7b_gates.csv
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R"))
X <- fread(file.path(TABDIR, "V4_S7b_dp_rates.csv"))[dgp == "GLMM-native" & arm == "H0" & method %in% MARGINAL_METHODS]
T1 <- fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))[dgp == "GLMM-native" & label == "manual", .(dgp, scen, ntime, delta = truth_marginal)]
T2 <- fread(file.path(TABDIR, "V4_S6b_ns_truth.csv"))[dgp == "GLMM-native" & label == "manual", .(dgp, scen, ntime, delta = tm)]
TT <- rbind(T1, T2)
X <- merge(X, TT, by = c("dgp", "scen", "ntime"), all.x = TRUE)
stopifnot(!anyNA(X$delta), nrow(X) == (480L + 24L + 192L) * 6L)
X[, np := NP[scen]]
exc <- function(d, s, cc) 100 * (pnorm(-cc + d / s) + pnorm(-cc - d / s) - 2 * pnorm(-cc))
X[, `:=`(ratio = delta / empse, cn = qnorm(0.975), ct = qt(0.975, df = N - np))]
X[, `:=`(excess_n = exc(delta, empse, cn), excess_t = exc(delta, empse, ct))]
X[n < 2, `:=`(excess_n = NA_real_, excess_t = NA_real_)]
R <- fread(file.path(TABDIR, "V4_S7b_rule_dp.csv"))[dgp == "GLMM-native" & method %in% MARGINAL_METHODS & ref %in% c("n", "t")]
R <- merge(R, X[, .(dgp, scen, N, ntime, method, delta, empse, ratio, excess_n, excess_t)], by = c("dgp", "scen", "N", "ntime", "method"))
R[, excess := fifelse(ref == "n", excess_n, excess_t)]
R[, rate_adj := rate - excess]
R[, hi_adj := assess & rate_adj - 2 * mc > 7.5]
R[, lo_adj := assess & rate_adj + 2 * mc < 2.5]
R[, group := fifelse(scen %in% c("M6", "M7"), "M6-M7", fifelse(scen == "M1c", "M1c", "M1-M5"))]
fwrite(R[, .(dgp, scen, group, N, ntime, src, method, ref, strategy, n, nd, rate, mc, delta, empse, ratio, excess, rate_adj, assess, hi, hi_adj, lo, lo_adj)],
       file.path(TABDIR, "V4_S7b_null_shift.csv"))
SM <- rbind(R[, .(n_dp = .N, max_excess = max(excess, na.rm = TRUE), mean_excess = mean(excess, na.rm = TRUE), max_ratio = max(abs(ratio), na.rm = TRUE),
                 n_hi = sum(hi), n_hi_adj = sum(hi_adj), n_hi_flip = sum(hi != hi_adj), n_lo_flip = sum(lo != lo_adj)), by = .(group, ref, method)],
            R[, .(n_dp = .N, max_excess = max(excess, na.rm = TRUE), mean_excess = mean(excess, na.rm = TRUE), max_ratio = max(abs(ratio), na.rm = TRUE),
                 n_hi = sum(hi), n_hi_adj = sum(hi_adj), n_hi_flip = sum(hi != hi_adj), n_lo_flip = sum(lo != lo_adj)), by = .(ref, method)][, group := "all"],
            R[, .(n_dp = .N, max_excess = max(excess, na.rm = TRUE), mean_excess = mean(excess, na.rm = TRUE), max_ratio = max(abs(ratio), na.rm = TRUE),
                 n_hi = sum(hi), n_hi_adj = sum(hi_adj), n_hi_flip = sum(hi != hi_adj), n_lo_flip = sum(lo != lo_adj)), by = .(group, ref)][, method := "all"], fill = TRUE)
setorder(SM, ref, group, method); fwrite(SM, file.path(TABDIR, "V4_S7b_null_summary.csv")); print(SM, digits = 3)
FL <- R[hi != hi_adj | lo != lo_adj]; fwrite(FL, file.path(TABDIR, "V4_S7b_null_flips.csv")); cat("flags that change after subtracting the excess:", nrow(FL), "\n")
cat("largest excess (normal):\n"); print(R[ref == "n"][order(-excess)][1:5, .(scen, N, ntime, method, delta, empse, rate, excess)], digits = 3)
# ---- GC4: Monte Carlo counter-check on 12 cells (4 largest excess + 8 at random, seed 20261021) --------------------
Xn <- X[is.finite(excess_n)]
set.seed(20261021)
pick <- unique(c(order(-Xn$excess_n)[1:4], sample(setdiff(seq_len(nrow(Xn)), order(-Xn$excess_n)[1:4]), 8)))
MC <- Xn[pick, .(scen, N, ntime, method, delta, empse, ratio, cn, ct, excess_n, excess_t)]
MC[, `:=`(mc_n = NA_real_, mc_t = NA_real_)]
for (i in seq_len(nrow(MC))) {
  z <- rnorm(1e6) + MC$ratio[i]
  MC$mc_n[i] <- 100 * (mean(abs(z) > MC$cn[i]) - 2 * pnorm(-MC$cn[i]))
  MC$mc_t[i] <- 100 * (mean(abs(z) > MC$ct[i]) - 2 * pnorm(-MC$ct[i]))
}
MC[, `:=`(dev_n = abs(mc_n - excess_n), dev_t = abs(mc_t - excess_t))]
gc4 <- MC[, all(dev_n < 0.05 & dev_t < 0.05)]
fwrite(MC, file.path(TABDIR, "V4_S7b_null_gc4.csv")); print(MC, digits = 3)
GT <- fread(file.path(TABDIR, "V4_S7b_gates.csv"))[gate != "GC4"]
GT <- rbind(GT, data.table(gate = "GC4", pass = gc4, detail = sprintf("12 cells (4 with the largest excess, 8 at random, seed 20261021), 1e6 normal draws each: max deviation from the formula %.3f points (normal), %.3f (t)", max(MC$dev_n), max(MC$dev_t))))
setorder(GT, gate); fwrite(GT, file.path(TABDIR, "V4_S7b_gates.csv")); print(GT)
