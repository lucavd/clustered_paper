# =============================================================================
# s7d2_explain.R -- explanations of the S7d2 results that differ from the predictions (written after the run):
# (1) finite-interval share (PG5, PG6) by T and N, H1, and why intervals are not finite: no estimate (|estimate| >= 10
#     or diverged full fit), unbounded with non-computable values of b0 met on the way (D-57 rules), unbounded without;
# (2) width ratio (PG7): interval of the test / GEE independence + Fay-Graubard + t, and GEE independence / GEE AR(1)
#     (both fewclust FG, t at N - p) from the saved S7d / S7e rows, H1, same datasets;
# (3) gate GG4, the intervals whose point 1e-4 beyond a finite end was not computable: reason at that point and p at
#     the bracketing rejected point.
# Outputs: tables/V4_S7d2_explain_finite.csv, V4_S7d2_explain_width.csv, V4_S7d2_GG4_detail.csv
# =============================================================================
source(path.expand("~/MMM/R/s7d2/s7d2_common.R")); setDTthreads(8L)
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
NPc <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
TK <- s7d2_tasks()[label == "mmlb"]
rd <- function(r) { t <- TK[r]; x <- readRDS(t$file); o <- readRDS(t$old_file); g <- if (t$dgp == DGP_COP) o$gee else o$rows
  merge(x, g[cor == "ind", .(sim_id, ind_ok = gee_ok, ind_est = est, ind_se_fg = se_fg)], by = "sim_id", all.x = TRUE) }
X <- rbindlist(parallel::mclapply(seq_len(nrow(TK)), rd, mc.cores = 16L), fill = TRUE)
X[, nd := !design_degenerate]; X[, np := NPc[scen]]; X[, qg := qt(0.975, N - np)]
X[, ui := nd & is.finite(est) & is.na(est_why)]; s9_drop_const(X, "ui", what = "interval of the test"); X[, fin := ui & is.finite(lo) & is.finite(hi)]
X[, ug := nd & ar1_gee_ok %in% TRUE & is.finite(ar1_est) & abs(ar1_est) < 10 & is.finite(ar1_se_fg) & ar1_se_fg > 0 & ar1_se_fg < 10]
X[, uind := nd & ind_ok %in% TRUE & is.finite(ind_est) & abs(ind_est) < 10 & is.finite(ind_se_fg) & ind_se_fg > 0 & ind_se_fg < 10]
s9_drop_const(X, "ug", what = "GEE AR(1) FG"); s9_drop_const(X, "uind", what = "GEE independence FG")   # D-82
X[, Ng := cut(N, c(0, 10, 25, 65, 145), labels = c("10", "15-25", "35-65", "75-145"))]
# (1) finite share and reasons, H1
F1 <- X[nd == TRUE, .(datasets = .N, finite = 100 * mean(fin), no_estimate = 100 * mean(!ui),
                      unbounded_with_noncomputable = 100 * mean(ui & !fin & n_na > 0), unbounded_without = 100 * mean(ui & !fin & n_na == 0),
                      gee_ar1_fg_usable = 100 * mean(ug)), by = .(ntime, Ng)][order(ntime, Ng)]
print(F1, digits = 3); fwrite(F1, file.path(TABDIR, "V4_S7d2_explain_finite.csv"))
# (2) width ratios on datasets where all three intervals are finite (test interval of width > 0), per point geometric mean, then mean of logs over points
X[, w_f := fifelse(fin, hi - lo, NA_real_)]; X[, w_i := fifelse(uind, 2 * qg * ind_se_fg, NA_real_)]; X[, w_a := fifelse(ug, 2 * qg * ar1_se_fg, NA_real_)]
P <- X[fin & uind & ug & w_f > 0, .(n = .N, l_fi = mean(log(w_f / w_i)), l_ia = mean(log(w_i / w_a)), l_fa = mean(log(w_f / w_a))), by = .(dgp, scen, N, ntime, Ng)][n >= 20]
W <- rbind(P[, .(group = "all points", points = .N, flip_over_ind = exp(mean(l_fi)), ind_over_ar1 = exp(mean(l_ia)), flip_over_ar1 = exp(mean(l_fa)))],
           P[, .(group = paste0("T = ", ntime), points = .N, flip_over_ind = exp(mean(l_fi)), ind_over_ar1 = exp(mean(l_ia)), flip_over_ar1 = exp(mean(l_fa))), by = ntime][order(ntime)][, -"ntime"],
           P[, .(group = paste0("N ", Ng), points = .N, flip_over_ind = exp(mean(l_fi)), ind_over_ar1 = exp(mean(l_ia)), flip_over_ar1 = exp(mean(l_fa))), by = Ng][order(Ng)][, -"Ng"],
           P[, .(group = dgp, points = .N, flip_over_ind = exp(mean(l_fi)), ind_over_ar1 = exp(mean(l_ia)), flip_over_ar1 = exp(mean(l_fa))), by = dgp][, -"dgp"])
print(W, digits = 3); fwrite(W, file.path(TABDIR, "V4_S7d2_explain_width.csv"))
# (3) GG4 cases (both arms)
TA <- s7d2_tasks()
Y <- rbindlist(parallel::mclapply(TA$file, function(f) { x <- readRDS(f)
  x[is.finite(est) & is.na(est_why) & ((is.finite(lo) & is.na(p_lo_out)) | (is.finite(hi) & is.na(p_hi_out)))] }, mc.cores = 16L), fill = TRUE)
G4 <- rbindlist(lapply(seq_len(nrow(Y)), function(r) { b <- Y[r]
  t <- unique(TA[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen); set.seed(FLIP_SEED(b$sim_id)); S <- flip_signs(P$N)
  side <- if (is.finite(b$lo) && is.na(b$p_lo_out)) "lower" else "upper"
  e <- if (side == "lower") b$lo else b$hi; rj <- if (side == "lower") b$rej_lo else b$rej_hi; out <- if (side == "lower") e - CI_EPS_BEYOND else e + CI_EPS_BEYOND
  data.table(dgp = b$dgp, scen = b$scen, N = b$N, ntime = b$ntime, label = b$label, sim_id = b$sim_id, side = side, end = e, rej = rj,
             p_end = flip_p_b(P, S, e)$p, p_rej = flip_p_b(P, S, rj)$p, why_beyond = flip_p_b(P, S, out)$why) }))
print(G4[, .N, by = .(why_beyond, rej_rejected = p_rej < CI_ALPHA, end_accepted = is.na(p_end) | p_end >= CI_ALPHA)])
print(G4[, .N, by = .(N, ntime)][order(N, ntime)])
fwrite(G4, file.path(TABDIR, "V4_S7d2_GG4_detail.csv"))
