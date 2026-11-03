# =============================================================================
# s7d2_agg.R -- S7d2 aggregation (written and committed before the production run): gates GG1, GG2, GG4; per point x
# arm coverage of the test at the true value and of the reported interval, share with a finite interval, width, width
# ratio with the GEE AR(1) + fewclust Fay-Graubard + t (N - p), coverage of that GEE; MCSE and Bradley-type flags;
# block-bootstrap aggregates (B = 2000, indices of results_v4/S4/boot_idx.rds and results_v4/S6b/ns/boot_idx_67.rds);
# predictions PG1-PG8 with their verdicts.
# Inputs : results_v4/S7d2/rows/*.rds; results_v4/S7d2/md5_results_full_{before,after}.txt
# Outputs: tables/V4_S7d2_gates.csv, V4_S7d2_dp.csv, V4_S7d2_boot_means.csv, V4_S7d2_predictions.csv,
#          V4_S7d2_by_T.csv, V4_S7d2_small.csv, V4_S7d2_reasons.csv; results_v4/S7d2/agg/boot_mats.rds
# With S7D2_DRY=1 (dry run) the tables go to results_v4/S7d2/agg/dry/ instead of tables/.
# =============================================================================
source(path.expand("~/MMM/R/s7d2/s7d2_common.R")); setDTthreads(8L)
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
t0 <- Sys.time(); dir.create(file.path(OUT_S7D2, "agg"), showWarnings = FALSE)
if (nzchar(Sys.getenv("S7D2_DRY"))) { TABDIR <- file.path(OUT_S7D2, "agg/dry"); dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE) }
NPc <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
TK <- s7d2_tasks(); have <- file.exists(TK$file)
X <- rbindlist(parallel::mclapply(TK$file[have], readRDS, mc.cores = 16L), fill = TRUE)
cat(sprintf("read %d chunks, %d replicates (%.1f min)\n", sum(have), nrow(X), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
# ---- GG1, GG2, GG4 ---------------------------------------------------------------------------------------------------
md5 <- function(tag) { f <- file.path(OUT_S7D2, sprintf("md5_results_full_%s.txt", tag)); if (file.exists(f)) readLines(f) else NA_character_ }
m_b <- md5("before"); m_a <- md5("after"); m_s7b <- readLines(file.path(REPO, "results_v4/S7b/md5_results_full_after.txt"))
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
Gn <- X[generated == TRUE]
E <- X[generated == TRUE & is.finite(est) & is.na(est_why)]
gg4 <- E[, .(n = .N, bad_order = sum(!(lo <= est & est <= hi)),
             bad_end_in = sum((is.finite(lo) & !is.na(p_lo) & p_lo < CI_ALPHA) | (is.finite(hi) & !is.na(p_hi) & p_hi < CI_ALPHA)),
             bad_end_out = sum((is.finite(lo) & !(!is.na(p_lo_out) & p_lo_out < CI_ALPHA)) | (is.finite(hi) & !(!is.na(p_hi_out) & p_hi_out < CI_ALPHA))),
             na_end_out = sum((is.finite(lo) & is.na(p_lo_out)) | (is.finite(hi) & is.na(p_hi_out))))]
GT <- data.table(
  gate = c("GG1", "GG1", "GG1", "GG1", "GG1", "GG2", "GG2", "GG2", "GG2", "GG4", "GG4", "GG4", "GG4"),
  check = c("chunk files", "replicates", "unique replicate keys", "md5 before = after", "md5 before = S7b list",
            "n_obs and sum_y equal to S7d/S7e", "p0 s7d_flip() = saved p_flip", "p0 new function = saved p_flip", "CI errors",
            "intervals checked", "lo <= est <= hi violated", "p < 0.05 at a reported end", "not rejected 1e-4 beyond a finite end"),
  value = c(sum(have), nrow(X), nrow(unique(X[, ..K])), identical(m_b, m_a), identical(sort(m_b), sort(m_s7b)),
            Gn[, sum(!(n_obs == old_n_obs & sum_y == old_sum_y), na.rm = TRUE) + sum(is.na(old_n_obs))],
            Gn[, sum(!same(p0_old, old_p_flip))], Gn[, sum(!same(p0_new, old_p_flip))], X[, sum(!is.na(ci_error))],
            gg4$n, gg4$bad_order, gg4$bad_end_in, gg4$bad_end_out),
  target = c(47520, 4752000, 4752000, 1, 1, 0, 0, 0, 0, NA, 0, 0, 0))
GT[, pass := is.na(target) | value == target]
GT <- rbind(GT, data.table(gate = "GG4", check = c("  of which p not computable 1e-4 beyond", "reported: finite intervals of width 0"),
                           value = c(gg4$na_end_out, E[, sum(is.finite(lo) & is.finite(hi) & hi == lo)]), target = NA, pass = NA))
print(GT); fwrite(GT, file.path(TABDIR, "V4_S7d2_gates.csv"))
# ---- per-replicate indicators ----------------------------------------------------------------------------------------
X[, nd := !design_degenerate]
X[, np := NPc[scen]]
X[, ut := nd & is.finite(p_theta)]; s9_drop_const(X, "ut", what = "test at the truth"); X[, cov_t := ut & p_theta >= CI_ALPHA]   # D-82
X[, ui := nd & is.finite(est) & is.na(est_why)]; s9_drop_const(X, "ui", what = "interval of the test"); X[, cov_i := ui & lo <= theta & theta <= hi]
X[, fin := ui & is.finite(lo) & is.finite(hi)]; X[, wid := fifelse(fin, hi - lo, NA_real_)]
X[, ug := nd & ar1_gee_ok %in% TRUE & is.finite(ar1_est) & abs(ar1_est) < 10 & is.finite(ar1_se_fg) & ar1_se_fg > 0 & ar1_se_fg < 10]
s9_drop_const(X, "ug", what = "GEE AR(1) FG")   # D-82
X[, qg := qt(0.975, N - np)]
X[, cov_g := ug & abs(ar1_est - theta) <= qg * ar1_se_fg]; X[, wid_g := fifelse(ug, 2 * qg * ar1_se_fg, NA_real_)]
# (after the run) 51 finite intervals have width 0 (centre not computable, neighbours rejected, p constant in b0):
# they stay in the finite share but not in the width ratio, where log(0) would enter as -Inf (s7d2_zerowidth.R)
X[, zw := fin & wid == 0]
X[, pr := fin & ug & wid > 0]; X[, lr := fifelse(pr, log(wid / wid_g), NA_real_)]
X[, small := N <= 25 & ntime <= 3]
# test at b0 = 0 as in S7d (saved p_flip) vs the same test with the D-57 rules (score numerically zero / diverged fit)
X[, u_old := nd & is.finite(old_p_flip)]; X[, u_g0 := nd & is.finite(p0_guard)]; X[, noise0 := u_old & !u_g0]
# ---- per point x arm -------------------------------------------------------------------------------------------------
pm <- function(k, n) 100 * k / n; ms <- function(p, n) sqrt(p * (100 - p) / n)
D <- X[, .(n_gen = sum(generated), nd = sum(nd), n_ut = sum(ut), cov_t = pm(sum(cov_t), sum(ut)), n_ui = sum(ui), cov_i = pm(sum(cov_i), sum(ui)),
           n_fin = sum(fin), pfin = pm(sum(fin), sum(nd)), wid_mean = mean(wid, na.rm = TRUE), wid_med = median(wid, na.rm = TRUE),
           n_ug = sum(ug), cov_g = pm(sum(cov_g), sum(ug)), pu_g = pm(sum(ug), sum(nd)), wid_g_med = median(wid_g, na.rm = TRUE),
           n_pr = sum(pr), mlr = mean(lr, na.rm = TRUE), sdlr = sd(lr, na.rm = TRUE), n_eval = mean(n_eval, na.rm = TRUE), n_na = sum(n_na, na.rm = TRUE),
           secs = mean(secs, na.rm = TRUE), pu_flip_s7d = pm(sum(u_old), sum(nd)), pu_flip_guard = pm(sum(u_g0), sum(nd)),
           rate_flip_guard = pm(sum(u_g0 & p0_guard < CI_ALPHA), sum(u_g0)), n_noise0 = sum(noise0)), by = .(dgp, scen, N, ntime, label, small)]
D[, `:=`(mcse_cov_t = ms(cov_t, n_ut), mcse_cov_i = ms(cov_i, n_ui), mcse_pfin = ms(pfin, nd), mcse_cov_g = ms(cov_g, n_ug),
         gm_ratio = exp(mlr), mcse_mlr = sdlr / sqrt(n_pr))]
D[, `:=`(lo_t = n_ut >= 20 & cov_t + 2 * mcse_cov_t < 92.5, hi_t = n_ut >= 20 & cov_t - 2 * mcse_cov_t > 97.5,
         lo_g = n_ug >= 20 & cov_g + 2 * mcse_cov_g < 92.5, hi_g = n_ug >= 20 & cov_g - 2 * mcse_cov_g > 97.5)]
DRY <- nzchar(Sys.getenv("S7D2_DRY"))   # dry run: partial grid, skip the completeness assertion
if (!DRY) stopifnot(nrow(D) == 4752L, D[label == "mmlb" & small == TRUE, .N] == 216L)
fwrite(D, file.path(TABDIR, "V4_S7d2_dp.csv"))
BT <- D[, .(points = .N, cov_t = mean(cov_t, na.rm = TRUE), cov_i = mean(cov_i, na.rm = TRUE), pfin = mean(pfin), wid_med = median(wid_med, na.rm = TRUE),
            gm_ratio = exp(mean(mlr, na.rm = TRUE)), cov_g = mean(cov_g, na.rm = TRUE), pu_g = mean(pu_g), n_lo_t = sum(lo_t), n_hi_t = sum(hi_t),
            pu_flip_s7d = mean(pu_flip_s7d), pu_flip_guard = mean(pu_flip_guard), rate_flip_guard = mean(rate_flip_guard, na.rm = TRUE)),
        by = .(mech = fifelse(dgp == DGP_COP, "copula", "study mechanisms"), label, ntime)][order(mech, label, ntime)]
fwrite(BT, file.path(TABDIR, "V4_S7d2_by_T.csv"))
SM <- D[, .(points = .N, pu_flip_s7d = mean(pu_flip_s7d), pu_flip_guard = mean(pu_flip_guard), rate_flip_guard = mean(rate_flip_guard, na.rm = TRUE),
            pfin = mean(pfin), pu_g = mean(pu_g), cov_t = mean(cov_t, na.rm = TRUE), n_noise0 = sum(n_noise0)), by = .(label, small)][order(label, small)]
SM[, set := fifelse(small, "small designs (T 2-3, N <= 25)", "other designs")]
SM <- rbind(SM, D[, .(small = NA, points = .N, pu_flip_s7d = mean(pu_flip_s7d), pu_flip_guard = mean(pu_flip_guard), rate_flip_guard = mean(rate_flip_guard, na.rm = TRUE),
            pfin = mean(pfin), pu_g = mean(pu_g), cov_t = mean(cov_t, na.rm = TRUE), n_noise0 = sum(n_noise0), set = "all 2,376 points"), by = label])
print(SM, digits = 4); fwrite(SM, file.path(TABDIR, "V4_S7d2_small.csv"))
RS <- rbind(X[nd == TRUE & !ui, .(what = "no interval", reason = fifelse(is.na(est_why), "not computed", est_why))][, .N, by = .(what, reason)],
            X[nd == TRUE & !ut, .(what = "test at truth not usable", reason = fifelse(is.na(theta_why), "NA", theta_why))][, .N, by = .(what, reason)],
            X[ui == TRUE & !fin, .(what = "interval not finite", reason = fifelse(is.finite(lo), "upper unbounded", fifelse(is.finite(hi), "lower unbounded", "both unbounded")))][, .N, by = .(what, reason)])
fwrite(RS, file.path(TABDIR, "V4_S7d2_reasons.csv")); print(RS)
# ---- block bootstrap -------------------------------------------------------------------------------------------------
BI <- c(readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds")), readRDS(file.path(REPO, "results_v4/S6b/ns/boot_idx_67.rds"))); names(BI) <- as.character(1:7)
CNT <- lapply(BI, function(Mx) { Cm <- matrix(0, 2000L, 1000L); for (b in 1:2000) Cm[b, ] <- tabulate(Mx[b, ], 1000L); t(rbind(1, Cm)) })
model_of <- function(scen) fifelse(scen == "M1c", 1L, as.integer(sub("^M", "", scen)))
DPK <- unique(X[, .(dgp, scen, N, ntime)])[order(dgp, scen, N, ntime)]; DPK[, key := paste(dgp, scen, N, ntime)]; DPK[, model := model_of(scen)]
X[, key := paste(dgp, scen, N, ntime)]
# ratio of weighted sums: 100 * sum(w num) / sum(w den), or mean of a value (scale = 1); 2,376 x 2001 (column 1 = estimate)
boot_ratio <- function(arm, num, den, scale = 100) {
  x0 <- X[label == arm]; out <- matrix(NA_real_, nrow(DPK), 2001L)
  for (m in unique(DPK$model)) {
    kk <- DPK[model == m]$key; x <- x0[key %in% kk]; ix <- cbind(match(x$key, kk), x$sim_id)
    Nu <- matrix(0, length(kk), 1000L); De <- Nu
    nv <- as.numeric(x[[num]]); nv[!is.finite(nv)] <- 0; Nu[ix] <- nv; De[ix] <- as.numeric(x[[den]])
    W <- CNT[[as.character(m)]]; out[match(kk, DPK$key), ] <- scale * (Nu %*% W) / (De %*% W)
  }
  out
}
summ <- function(M, sel = rep(TRUE, nrow(M)), f = identity) { mm <- f(colMeans(M[sel, , drop = FALSE], na.rm = TRUE))
  c(est = mm[1], lo = unname(quantile(mm[-1], 0.025, na.rm = TRUE)), hi = unname(quantile(mm[-1], 0.975, na.rm = TRUE)), sd = sd(mm[-1], na.rm = TRUE)) }
X[, lr0 := fifelse(pr, lr, 0)]
B <- list(cov_t_H1 = boot_ratio("mmlb", "cov_t", "ut"), cov_i_H1 = boot_ratio("mmlb", "cov_i", "ui"), pfin_H1 = boot_ratio("mmlb", "fin", "nd"),
          mlr_H1 = boot_ratio("mmlb", "lr0", "pr", scale = 1), cov_g_H1 = boot_ratio("mmlb", "cov_g", "ug"),
          cov_t_H0 = boot_ratio("manual", "cov_t", "ut"), pfin_H0 = boot_ratio("manual", "fin", "nd"), cov_g_H0 = boot_ratio("manual", "cov_g", "ug"))
sm <- DPK$N <= 25 & DPK$ntime <= 3; cop <- DPK$dgp == DGP_COP
BM <- rbindlist(list(
  data.table(q = "PG1 coverage of the test at the truth, H1", t(summ(B$cov_t_H1))),
  data.table(q = "PG4 coverage interval - coverage test, H1", t(summ(B$cov_i_H1 - B$cov_t_H1))),
  data.table(q = "PG5 finite interval share, H1", t(summ(B$pfin_H1))),
  data.table(q = "PG6 finite interval share, H1, small points", t(summ(B$pfin_H1, sm))),
  data.table(q = "PG7 width ratio flip / GEE AR(1) FG t (geometric), H1", t(summ(B$mlr_H1, f = exp))),
  data.table(q = "PG8 coverage GEE AR(1) FG t, H1", t(summ(B$cov_g_H1))),
  data.table(q = "coverage of the test at the truth, H0", t(summ(B$cov_t_H0))),
  data.table(q = "coverage GEE AR(1) FG t, H0", t(summ(B$cov_g_H0))),
  data.table(q = "finite interval share, H0", t(summ(B$pfin_H0))),
  data.table(q = "coverage of the test, H1, copula", t(summ(B$cov_t_H1, cop))),
  data.table(q = "coverage of the test, H1, study mechanisms", t(summ(B$cov_t_H1, !cop))),
  data.table(q = "width ratio, H1, small points", t(summ(B$mlr_H1, sm, f = exp))),
  data.table(q = "coverage test - coverage GEE AR(1) FG t, H1", t(summ(B$cov_t_H1 - B$cov_g_H1)))))
setnames(BM, c("q", "est", "lo", "hi", "sd")); print(BM, digits = 4); fwrite(BM, file.path(TABDIR, "V4_S7d2_boot_means.csv"))
saveRDS(list(DPK = DPK, B = B), file.path(if (DRY) TABDIR else file.path(OUT_S7D2, "agg"), "boot_mats.rds"))
# ---- predictions -----------------------------------------------------------------------------------------------------
verdict <- function(lo, hi, a, b) if (lo >= a && hi <= b) "confirmed" else if (hi < a || lo > b) "refuted" else "indeterminate"
# qq is deliberately not called q: inside BM[...] the name q would be the column, and every row would match
row_bm <- function(code, qq, a, b) { r <- BM[startsWith(BM$q, qq)]; stopifnot(nrow(r) == 1L); data.table(code = code, obs = r$est, lo = r$lo, hi = r$hi, a = a, b = b, verdict = verdict(r$lo, r$hi, a, b)) }
row_ct <- function(code, v, a, b) data.table(code = code, obs = v, lo = v, hi = v, a = a, b = b, verdict = verdict(v, v, a, b))
H1 <- D[label == "mmlb"]
PR <- rbindlist(list(row_bm("PG1", "PG1", 92.5, 97.5), row_ct("PG2", H1[, sum(lo_t)], 0, 0), row_ct("PG3", H1[, sum(hi_t)], 0, 60),
                     row_bm("PG4", "PG4", -1, 1), row_bm("PG5", "PG5", 90, 100), row_bm("PG6", "PG6", 60, 100),
                     row_bm("PG7", "PG7", 0.90, 1.10), row_bm("PG8", "PG8", 92.5, 97.5)))
print(PR, digits = 4); fwrite(PR, file.path(TABDIR, "V4_S7d2_predictions.csv"))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
