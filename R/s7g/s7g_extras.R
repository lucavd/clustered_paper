# =============================================================================
# s7g_extras.R -- S7g: numbers for the text points of the step (docs/steps/S7g.md, doubts 1, 2, 4, 6, 8, 9).
#   E1 (doubt 1) what "never clearly above 7.5%" guarantees: per Table 4 strategy, largest observed rate with its usable
#      replicates and MCSE, points clearly above, points where an excess is not excluded (exact upper limit > 7.5% but not
#      clearly above), points not assessed (< 20 usable); check: clearly-above counts equal tables/V5_S7f_table4.csv.
#      Fay-Graubard across working correlations: per-point rate differences (t reference), three study mechanisms.
#   E2 (doubt 2) denominators: the per-point denominators of the Table 3 chain (S7b) and of Table 4 (S7d) are compared.
#   E3 (doubt 4) pseudo-true value of the mixed model on Copula-AR1 (tables/V4_S7e_glmm_pseudotrue.csv), read back.
#   E4 (doubt 6) M6-M7, marginalized model with the check: bias / empirical SE, model SE / empirical SE, coverage (t), by
#      arm, and the same pooled separately for fits where the check intervened or not (description, not a criterion).
#   E5 (doubt 2, cohort) usable share over all 500 subsamples at the smallest N (tables/V5_S7f_lus_tab316.csv).
#   E6 (doubt 9) volume of computation: datasets and production-run wall time from the run logs of results_v4.
# Outputs: tables/V5_S7g_calib.csv, V5_S7g_fg_wc.csv, V5_S7g_m67_cc.csv, V5_S7g_runs.csv, V5_S7g_extra.csv (key, value, mcse,
#          kind, source, note: the registry source of the new tokens f:g7|<key>). Read-only on every input.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L); setwd("/home/user/MMM")
cp_lo <- function(k, n) fifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
cp_up <- function(k, n) fifelse(k >= n, 1, qbeta(0.975, k + 1, pmax(n - k, 1e-12)))
EX <- list(); put <- function(key, value, kind = "sig2", source = "", note = "", mcse = NA_real_)
  { x <- data.table(k_ = key, value = as.character(value), mcse = if (is.na(mcse)) "" else sprintf("%.17g", mcse), kind = kind, source = source, note = note)
    setnames(x, "k_", "key"); EX[[length(EX) + 1]] <<- x }
## ---- E1 ----------------------------------------------------------------------------------------------------------------
A  <- fread("tables/V5_S7f_cp_S7d_dp_rates.csv"); B <- fread("tables/V5_S7f_cp_S7e_dp_rates.csv")
C  <- fread("tables/V5_S7f_cp_S7b_rule_dp.csv"); FL <- fread("tables/V5_S7f_cp_S7d2_flip_rules_dp.csv")
est <- c("gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_ar1_md", "gee_ind_md", "gee_exch_md", "gee_ar1_kc", "gee_ind_kc",
         "gee_exch_kc", "gee_ind_cr2", "gee_exch_cr2", "gee_ar1_unc", "gee_ind_unc", "gee_exch_unc")
cols <- c("dgp", "scen", "N", "ntime", "strategy", "n", "rate", "mc", "assess", "hi", "lo")
s1 <- A[ref == "t" & strategy %in% c(est, "study_mmm"), ..cols][strategy == "study_mmm", strategy := "mmm_cc"]
s2 <- C[(ref == "t" & method %in% c("glmer", "shrink")) | (ref == "n" & method == "glmer_lrt"), .(dgp, scen, N, ntime, strategy = method, n, rate, mc, assess, hi, lo)]
s3 <- B[(ref == "t" & strategy %in% c(est, "mmm_cc", "glmer", "shrink")) | (ref == "n" & strategy == "glmer_lrt"), ..cols]
s4 <- FL[version == "with D-57 rules" & label == "manual", .(dgp, scen, N, ntime, strategy = "flip", n, rate, mc = mcse, assess = n >= 20, hi, lo)]
D <- rbind(s1, s2, s3, s4, use.names = TRUE)
stopifnot(D[, .N, by = strategy][, all(N == 2376L)])
D[, k := fifelse(n > 0, round(n * rate / 100), 0)]
stopifnot(D[n > 0, max(abs(k - n * rate / 100))] < 1e-6)
D[, up := 100 * cp_up(k, n)]
D[, notexcl := assess & !hi & up > 7.5]
T4 <- fread("tables/V5_S7f_table4.csv")
CAL <- D[, { a <- .SD[assess == TRUE]; j <- which.max(a$rate)
  .(n_points = .N, n_hi = sum(hi), n_notexcl = sum(notexcl), n_notassessed = sum(!assess), max_rate = a$rate[j], max_n = a$n[j],
    max_mc = a$mc[j], max_where = sprintf("%s %s N=%d T=%d", a$dgp[j], a$scen[j], a$N[j], a$ntime[j]), median_n = as.numeric(median(n)),
    max_up_not_hi = if (any(!a$hi)) max(a$up[!a$hi]) else NA_real_) }, by = strategy]
CAL <- merge(CAL, T4[, .(strategy, hi_all)], by = "strategy")
stopifnot(all(CAL$n_hi == CAL$hi_all))
fwrite(CAL, "tables/V5_S7g_calib.csv"); print(CAL[order(n_hi)], digits = 3)
for (s in CAL$strategy) { r <- CAL[strategy == s]
  put(paste0("cal_maxrate_", s), r$max_rate, "sig2", "tables/V5_S7g_calib.csv", "largest rate among assessed points", r$max_mc)
  put(paste0("cal_maxn_", s), r$max_n, "int", "tables/V5_S7g_calib.csv", "usable replicates at that point")
  put(paste0("cal_notexcl_", s), r$n_notexcl, "int", "tables/V5_S7g_calib.csv", "assessed, not clearly above, exact upper limit > 7.5%")
  put(paste0("cal_notass_", s), r$n_notassessed, "int", "tables/V5_S7g_calib.csv", "points with < 20 usable") }
thr <- function(n) { k <- 0:n; min(k[100 * cp_lo(k, n) > 7.5]) }
put("cal_thr1000", 100 * thr(1000) / 1000, "plain", "computed", "smallest rate flagged clearly above with 1000 usable")
put("cal_thr500", 100 * thr(500) / 500, "plain", "computed", "smallest rate flagged clearly above with 500 usable")
put("cal_thr100", 100 * thr(100) / 100, "plain", "computed", "smallest rate flagged clearly above with 100 usable")
## Fay-Graubard across working correlations (three study mechanisms, t reference, points assessed for all three)
W <- dcast(A[ref == "t" & strategy %in% c("gee_exch_fg", "gee_ind_fg", "gee_ar1_fg")], dgp + scen + N + ntime ~ strategy, value.var = c("rate", "assess"))
W <- W[assess_gee_exch_fg & assess_gee_ind_fg & assess_gee_ar1_fg]
W[, `:=`(d_ind = rate_gee_ind_fg - rate_gee_exch_fg, d_ar1 = rate_gee_ar1_fg - rate_gee_exch_fg)]
FG <- rbindlist(lapply(c("d_ind", "d_ar1"), function(v) { x <- W[[v]]
  data.table(contrast = v, n_points = length(x), mean = mean(x), q05 = quantile(x, 0.05), median = median(x), q95 = quantile(x, 0.95),
             max_abs = max(abs(x)), n_abs_gt_2.5 = sum(abs(x) > 2.5), where_max = W[which.max(abs(x)), sprintf("%s %s N=%d T=%d", dgp, scen, N, ntime)]) }))
fwrite(FG, "tables/V5_S7g_fg_wc.csv"); print(FG, digits = 3)
for (v in c("d_ind", "d_ar1")) { r <- FG[contrast == v]
  put(paste0("fgwc_npts"), r$n_points, "int", "tables/V5_S7g_fg_wc.csv", "points assessed for all three working correlations")
  put(paste0("fgwc_q05_", v), r$q05, "sig2", "tables/V5_S7g_fg_wc.csv"); put(paste0("fgwc_q95_", v), r$q95, "sig2", "tables/V5_S7g_fg_wc.csv")
  put(paste0("fgwc_max_", v), r$max_abs, "sig2", "tables/V5_S7g_fg_wc.csv"); put(paste0("fgwc_gt25_", v), r$n_abs_gt_2.5, "int", "tables/V5_S7g_fg_wc.csv") }
EX <- EX[!duplicated(vapply(EX, `[[`, "", "key"), fromLast = TRUE)]
## ---- E2 ----------------------------------------------------------------------------------------------------------------
dn <- merge(A[ref == "t" & strategy == "gee_ar1_fg", .(dgp, scen, N, ntime, nd_t4 = nd)], C[strategy == "glmer|t", .(dgp, scen, N, ntime, nd_s7b = nd)], by = c("dgp", "scen", "N", "ntime"))
DEG <- fread("tables/V4_S2_degenerate_designs.csv")
put("den_points_same", dn[, sum(nd_t4 == nd_s7b)], "int", "V5_S7f_cp_S7d_dp_rates.csv vs V5_S7f_cp_S7b_rule_dp.csv", sprintf("of %d points", nrow(dn)))
put("den_deg_total", nrow(DEG), "int", "tables/V4_S2_degenerate_designs.csv", "degenerate (model, N, T, replicate) designs")
## ---- E3 ----------------------------------------------------------------------------------------------------------------
PT <- fread("tables/V4_S7e_glmm_pseudotrue.csv"); print(PT)
p1 <- PT[model == 1L]
put("pt_nsub", 20000L, "int", "R/s7e/s7e_glmm_pseudotrue.R", "subjects per dataset")
put("pt_min", min(p1$est_timetrt), "sig2", "tables/V4_S7e_glmm_pseudotrue.csv", "random intercept, T 2/4/10")
put("pt_max", max(p1$est_timetrt), "sig2", "tables/V4_S7e_glmm_pseudotrue.csv")
put("pt_semax", max(p1$se_timetrt), "sig2", "tables/V4_S7e_glmm_pseudotrue.csv")
## ---- E4 ----------------------------------------------------------------------------------------------------------------
V <- readRDS("results_v4/S6b/ns/rows_ns.rds")[method == "mmm"]
cat("rows_ns columns:", paste(names(V), collapse = ","), "\n")
if (!"th" %in% names(V)) { TR <- fread("tables/V4_S6b_ns_truth.csv"); V <- merge(V, TR[, .(dgp, scen, ntime, label, th = tm)], by = c("dgp", "scen", "ntime", "label")) }
if (!"deg" %in% names(V)) V[, deg := design_degenerate %in% TRUE]
if ("term" %in% names(V)) V <- V[!is.na(term)]
if (!"U0" %in% names(V)) stop("rows_ns.rds has no U0 column")
CC <- rbindlist(lapply(list.files("results_v4/S6b/ns/rows", "\\.rds$", full.names = TRUE), function(f) readRDS(f)$cc[, .(dgp, scen, N, ntime, label, sim_id, used, bnd)]))
V <- merge(V, CC, by = c("dgp", "scen", "N", "ntime", "label", "sim_id"), all.x = TRUE)
source("/home/user/MMM/R/s11/s11_cv.R"); s11_drop_cv(V, "U0", rows = rep(TRUE, nrow(V)), what = "rows_ns mmm")   # D-90
V[, u := U0 & !deg]; V[, arm := fifelse(label == "mmlb", "H1", "H0")]
V[, `:=`(err = estimate - th, cov_t = abs(estimate - th) <= qt(0.975, N - 4L) * se)]
PP <- V[u == TRUE, .(n = .N, bias = mean(err), empse = sd(estimate), modse = mean(se), kcov = sum(cov_t), share_B = mean(used == "B", na.rm = TRUE)),
        by = .(dgp, scen, N, ntime, arm)]
PP[, `:=`(rb = bias / empse, ser = modse / empse - 1, cov = 100 * kcov / n)]
PP[, `:=`(cov_hi = n >= 20 & 100 * cp_lo(n - kcov, n) > 7.5, cov_lo = n >= 20 & 100 * cp_up(n - kcov, n) < 2.5)]
S67 <- PP[n >= 20, .(n_points = .N, mean_rb = mean(rb), max_abs_rb = max(abs(rb)), median_ser = median(ser), min_ser = min(ser), max_ser = max(ser),
                     mean_cov = mean(cov), min_cov = min(cov), n_cov_clearly_low = sum(cov_hi), n_cov_clearly_high = sum(cov_lo), mean_share_B = 100 * mean(share_B)),
          by = .(dgp, scen, arm)]
# by intervention of the check: per design point (>= 20 usable fits in the group), then summarised over points (pooling
# estimates across N and T would inflate the empirical SE)
BYP <- V[u == TRUE & !is.na(used), .(n = .N, rb = mean(err) / sd(estimate), ser = mean(se) / sd(estimate) - 1, cov = 100 * mean(cov_t)),
         by = .(dgp, scen, N, ntime, arm, used)][n >= 20]
BYI <- BYP[, .(n_points = .N, fits = sum(n), mean_rb = mean(rb), max_abs_rb = max(abs(rb)), median_ser = median(ser), mean_cov = mean(cov)), by = .(dgp, scen, arm, used)]
BYI[, `:=`(rb = mean_rb, ser = median_ser, cov = mean_cov)]
fwrite(S67, "tables/V5_S7g_m67_cc.csv"); fwrite(BYI, "tables/V5_S7g_m67_cc_byused.csv"); fwrite(PP, "tables/V5_S7g_m67_cc_dp.csv")
print(S67, digits = 3); print(BYI, digits = 3)
h1 <- S67[arm == "H1"]
put("m67_rb_max", max(abs(S67$max_abs_rb)), "sig2", "tables/V5_S7g_m67_cc.csv", "largest |bias / empirical SE| over points, both arms")
put("m67_ser_min", min(S67$min_ser) * 100, "sig2", "tables/V5_S7g_m67_cc.csv", "smallest model SE / empirical SE - 1 (%) over points")
put("m67_ser_max", max(S67$max_ser) * 100, "sig2", "tables/V5_S7g_m67_cc.csv", "largest (%)")
put("m67_ser_med_min", min(S67$median_ser) * 100, "sig2", "tables/V5_S7g_m67_cc.csv", "smallest scenario median (%)")
put("m67_ser_med_max", max(S67$median_ser) * 100, "sig2", "tables/V5_S7g_m67_cc.csv", "largest scenario median (%)")
put("m67_cov_min", min(S67$mean_cov), "sig2", "tables/V5_S7g_m67_cc.csv", "smallest scenario mean coverage (%), t")
put("m67_cov_max", max(S67$mean_cov), "sig2", "tables/V5_S7g_m67_cc.csv", "largest scenario mean coverage (%), t")
put("m67_cov_low", sum(S67$n_cov_clearly_low), "int", "tables/V5_S7g_m67_cc.csv", "points with non-coverage clearly above 7.5%")
put("m67_cov_high", sum(S67$n_cov_clearly_high), "int", "tables/V5_S7g_m67_cc.csv", "points with non-coverage clearly below 2.5%")
put("m67_npts", S67[, sum(n_points)], "int", "tables/V5_S7g_m67_cc.csv", "point x arm assessed")
for (us in c("A", "B")) { b <- BYI[used == us]
  put(paste0("m67_cov_", us, "_min"), min(b$cov), "sig2", "tables/V5_S7g_m67_cc_byused.csv"); put(paste0("m67_cov_", us, "_max"), max(b$cov), "sig2", "tables/V5_S7g_m67_cc_byused.csv")
  put(paste0("m67_rb_", us, "_maxabs"), max(abs(b$rb)), "sig2", "tables/V5_S7g_m67_cc_byused.csv")
  put(paste0("m67_ser_", us, "_min"), 100 * min(b$ser), "sig2", "tables/V5_S7g_m67_cc_byused.csv"); put(paste0("m67_ser_", us, "_max"), 100 * max(b$ser), "sig2", "tables/V5_S7g_m67_cc_byused.csv") }
## ---- E5 ----------------------------------------------------------------------------------------------------------------
LT <- fread("tables/V5_S7f_lus_tab316.csv")[spec == "pieno" & N == min(N)]
for (m in c("GLMM", "MMM", "GEE", "GEE-FG", "GEE-MD")) put(paste0("lus_u500_", m), LT[method == m, pct_usable_500], "sig2", "tables/V5_S7f_lus_tab316.csv", "usable over all 500 subsamples, smallest N")
# draws without an invasively ventilated patient against the hypergeometric expectation (simple random draws, rev2/T2_defs.R)
LL <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds"); setDT(LL)
Kinv <- LL[, .(inv = any(invasive == 1)), by = id][, sum(inv)]; Ppat <- uniqueN(LL$id)
D3 <- fread("rev2/T3_degeneracy.csv")
DR <- D3[, .(draws = .N, none = sum(n_inv_pat == 0), mean_inv = mean(n_inv_pat)), by = N]
DR[, `:=`(p0 = dhyper(0, Kinv, Ppat - Kinv, N), exp_inv = N * Kinv / Ppat)]
DR[, z := fifelse(p0 > 0.001, (none - draws * p0) / sqrt(draws * p0 * (1 - p0)), NA_real_)]
DR[, p_two := mapply(function(x, n, p) binom.test(x, n, p)$p.value, none, draws, p0)]
fwrite(DR, "tables/V5_S7g_lus_draws.csv"); print(DR, digits = 3)
put("lus_draw_z10", DR[N == min(N), z], "sig2", "tables/V5_S7g_lus_draws.csv", "binomial z of draws with no ventilated patient, smallest N")
put("lus_draw_p10", DR[N == min(N), p_two], "sig2", "tables/V5_S7g_lus_draws.csv", "two-sided exact binomial p, smallest N")
put("lus_draw_zmax_other", DR[N != min(N) & is.finite(z), max(abs(z))], "sig2", "tables/V5_S7g_lus_draws.csv", "largest |z| at the other sizes with p0 > 0.001")
put("lus_draw_meaninv10", DR[N == min(N), mean_inv], "sig2", "tables/V5_S7g_lus_draws.csv", "mean ventilated patients per draw, smallest N")
put("lus_draw_expinv10", DR[N == min(N), exp_inv], "sig2", "tables/V5_S7g_lus_draws.csv", "expected under simple random sampling")
## ---- E6 ----------------------------------------------------------------------------------------------------------------
logs <- c(S3 = "results_v4/S3/run.log", S6b_ns = "results_v4/S6b/ns/run_progress.log", S7d = "results_v4/S7d/run_progress.log",
          S7e = "results_v4/S7e/run_progress.log", S7d2 = "results_v4/S7d2/run_progress.log", S7g = "results_v4/S7g/run_progress.log")
RU <- rbindlist(lapply(names(logs), function(k) { x <- readLines(logs[[k]]); d <- grep("DONE", x, value = TRUE)
  h <- vapply(d, function(z) { m <- regmatches(z, regexpr("[0-9.]+ (h|min)", z)); v <- as.numeric(sub(" .*", "", m)); if (grepl("min", m)) v / 60 else v }, 0)
  data.table(run = k, file = logs[[k]], done_lines = length(d), hours_last_segment = if (length(h)) tail(h, 1) else NA_real_, all_done_hours = paste(round(h, 2), collapse = ";")) }))
fwrite(RU, "tables/V5_S7g_runs.csv"); print(RU)
X <- rbindlist(EX); X <- X[!duplicated(key, fromLast = TRUE)]
fwrite(X, "tables/V5_S7g_extra.csv"); cat(sprintf("extra: %d keys\n", nrow(X)))
