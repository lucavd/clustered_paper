# =============================================================================
# s7e_agg.R -- S7e aggregation: gates GF1 (completeness, md5) and GF4 (b, c), per design point x arm x strategy rates with
# MCSE and Bradley flags on the copula mechanism, bias and coverage of the marginal estimators, the rule of Figure 1 on the
# 480 copula points and on the 2,376 points of all four mechanisms, block-bootstrap comparisons with the three study
# mechanisms at the same points (indices of results_v4/S4/boot_idx.rds), and predictions PF1-PF11 with their verdicts.
# Inputs : results_v4/S7e/rows/*.rds; results_v4/S7d/agg/boot_mats.rds (flip and GEE AR(1) FG of the study mechanisms);
#          results_v4/variant_cc/S4/rows_v4.rds (MMM with convergence check and GLMM Wald of the study mechanisms);
#          tables/V4_S7d_dp_rates.csv, tables/V4_S7b_rule_dp.csv (rule of Figure 1 on the 1,896 points).
# Outputs: tables/V4_S7e_GF1.csv, V4_S7e_GF4.csv, V4_S7e_dp_rates.csv, V4_S7e_bias_cov.csv, V4_S7e_rule.csv,
#          V4_S7e_boot_means.csv, V4_S7e_predictions.csv, V4_S7e_usable_reasons.csv
# =============================================================================
source(path.expand("~/MMM/R/s7e/s7e_common.R"))
source(path.expand("~/MMM/R/s7b/s7b_common.R"))       # TABDIR, NP, flags(), u0_flag(), add_usability()
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
source(path.expand("~/MMM/R/s11/s11_cv.R"))   # D-90 (S11, D-91): no result for the MMM with the check when its refit does not converge
AGG_CORES <- as.integer(Sys.getenv("S7E_AGG_CORES", "16")); data.table::setDTthreads(min(8L, AGG_CORES))
t0 <- Sys.time(); OUT <- OUT_S7E; dir.create(file.path(OUT, "agg"), showWarnings = FALSE)
TK <- s7e_tasks(); have <- file.exists(TK$file)
L <- parallel::mclapply(TK$file[have], readRDS, mc.cores = AGG_CORES)
RECS <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); ST <- rbindlist(lapply(L, `[[`, "study"), fill = TRUE)
CC <- rbindlist(lapply(L, `[[`, "cc"), fill = TRUE); GEE <- rbindlist(lapply(L, `[[`, "gee"), fill = TRUE)
SECS <- do.call(rbind, lapply(L, `[[`, "secs")); rm(L); invisible(gc())
cat(sprintf("read %d chunks: %d recs, %d study rows, %d cc, %d gee rows (%.1f min)\n", sum(have), nrow(RECS), nrow(ST), nrow(CC), nrow(GEE),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPC <- c("scen", "N", "ntime")
# ---- GF1 --------------------------------------------------------------------------------------------------------
gen_ok <- RECS[generated == TRUE, ..K6]
gee_n <- GEE[, .(cnt = .N), by = K6]
st_m <- merge(RECS[generated == TRUE, c(K6, "design_degenerate"), with = FALSE], ST[method %in% TECHS_S7E, .(nm = uniqueN(method)), by = K6], by = K6, all.x = TRUE)
st_m[is.na(nm), nm := 0L]   # a technique without a time:trt row and without a trace row returned a fit in which time:trt is not estimable
md5f <- file.path(OUT, "md5_results_full_after.txt"); md5b <- file.path(OUT, "md5_results_full_before.txt"); md5s <- file.path(REPO, "results_v4/S7b/md5_results_full_after.txt")
md5_ok <- file.exists(md5f) && identical(readLines(md5f), readLines(md5b)) && identical(readLines(md5f), readLines(md5s))
ge1 <- data.table(chunks = sum(have), chunks_expected = nrow(TK), recs = nrow(RECS), recs_unique = uniqueN(RECS[, ..K6]),
                  gen_errors = RECS[generated == FALSE, .N], gen_errors_nondegenerate = RECS[generated == FALSE & design_degenerate == FALSE, .N],
                  gee3_all = nrow(gee_n) == nrow(gen_ok) && all(gee_n$cnt == 3L), study4_all = all(st_m$nm == 4L), study4_missing = st_m[nm < 4L, .N], study4_missing_nondeg = st_m[nm < 4L & design_degenerate == FALSE, .N],
                  md5_unchanged = md5_ok)
ge1[, pass := chunks == chunks_expected & recs == 960000L & recs_unique == recs & gen_errors_nondegenerate == 0L & gee3_all & study4_all & md5_unchanged]
print(ge1); fwrite(ge1, file.path(TABDIR, "V4_S7e_GF1.csv"))
# ---- GF4 (b), (c) -------------------------------------------------------------------------------------------------
same <- function(a, b, rel = 0, abs_tol = 1e-10) abs(a - b) <= pmax(abs_tol, rel * pmax(abs(a), abs(b)))
b1 <- merge(ST[method == "gee" & term %in% "time:trt", c(K6, "estimate"), with = FALSE], GEE[cor == "exch", c(K6, "est"), with = FALSE], by = K6)
b1 <- b1[is.finite(estimate) & is.finite(est)]
b2 <- merge(ST[method == "mmm" & term %in% "time:trt", c(K6, "estimate"), with = FALSE], CC[A_ok == TRUE & A_has_term == TRUE, c(K6, "A_est"), with = FALSE], by = K6)
b2 <- b2[is.finite(estimate) & is.finite(A_est)]
cc4 <- GEE[cor == "ind" & sim_id <= 50 & gee_ok & is.finite(est) & abs(est) < 10]
c_md <- cc4[is.finite(se_md) & is.finite(smv_md) & se_md < 10]; c_fg <- cc4[is.finite(se_fg) & is.finite(smv_fg) & se_fg < 10]
GE4 <- data.table(part = c("b: study GEE exch = added GEE exch", "b: MMM fit A = study MMM (M3-M5)", "c: geesmv = fewclust MD (ind, 1-50)", "c: geesmv = fewclust FG (ind, 1-50)"),
                  n = c(nrow(b1), nrow(b2), nrow(c_md), nrow(c_fg)),
                  n_mismatch = c(b1[, sum(!same(estimate, est))], b2[, sum(!same(estimate, A_est))], c_md[, sum(!same(se_md, smv_md, 1e-6))], c_fg[, sum(!same(se_fg, smv_fg, 1e-6))]),
                  max_abs_diff = c(b1[, max(abs(estimate - est))], b2[, max(abs(estimate - A_est))], c_md[, max(abs(se_md - smv_md))], c_fg[, max(abs(se_fg - smv_fg))]))
GE4[, pass := n > 0 & n_mismatch == 0L]; print(GE4); fwrite(GE4, file.path(TABDIR, "V4_S7e_GF4.csv"))
# ---- per-replicate indicators (copula) -----------------------------------------------------------------------------
BASE <- RECS[, c(K6, "design_degenerate"), with = FALSE]; BASE[, np := NP[scen]]
truth <- function(scen, label) fifelse(label == "manual", 0, fifelse(scen == "M4", 0.05, 0.2))
S <- ST[method %in% c("mmm", "gee", "glmer", "shrink", "glmer_lrt") & term %in% "time:trt"]
Sm <- merge(S[method == "mmm"], CC[A_ok == TRUE & A_has_term == TRUE, c(K6, "cc_est", "cc_se", "cc_p", "cc_status", "used"), with = FALSE], by = K6, all.x = TRUE)
Mcc <- copy(Sm); Mcc[!is.na(cc_est), `:=`(estimate = cc_est, se = cc_se, p = cc_p, status = fifelse(used == "B", cc_status, status))]; Mcc[, method := "mmm_cc"]
Mst <- copy(Sm); Mst[, method := "mmm_study"]
S <- rbind(S[method != "mmm"], Mcc[, names(S), with = FALSE], Mst[, names(S), with = FALSE])
S[, U0 := u0_flag(method, converged, estimate, se, p)]
ind <- list()
for (m in unique(S$method)) {
  X <- merge(BASE, S[method == m, c(K6, "estimate", "se", "p", "U0"), with = FALSE], by = K6, all.x = TRUE)
  X[, use := !design_degenerate & U0 %in% TRUE]
  s9_drop_const(X, "use", what = m)   # D-82
  if (m == "mmm_cc") s11_drop_cv(X, "use", rows = rep(TRUE, nrow(X)), what = m)   # D-90
  if (m == "glmer_lrt") X[, `:=`(rej = use & p < 0.05, rejn = use & p < 0.05)]
  else X[, `:=`(rej = use & 2 * pt(-abs(estimate / se), N - np) < 0.05, rejn = use & 2 * pnorm(-abs(estimate / se)) < 0.05)]
  ind[[m]] <- X
}
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) {
  if (ty == "cr2" && cn == "ar1") next
  X <- merge(BASE, GEE[cor == cn, c(K6, "est", "gee_ok", SEC[[ty]]), with = FALSE], by = K6, all.x = TRUE); setnames(X, c("est", SEC[[ty]]), c("estimate", "se"))
  X[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(estimate) & abs(estimate) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(X, "use", what = paste0("gee_", cn, "_", ty))   # D-82
  X[, `:=`(rej = use & 2 * pt(-abs(estimate / se), N - np) < 0.05, rejn = use & 2 * pnorm(-abs(estimate / se)) < 0.05)]
  ind[[paste0("gee_", cn, "_", ty)]] <- X
}
for (v in c("flip", "flip_strict")) {
  X <- merge(BASE, RECS[, c(K6, if (v == "flip") "p_flip" else "p_flip_strict"), with = FALSE], by = K6); setnames(X, ncol(X), "p")
  X[, use := !design_degenerate & is.finite(p)]; s9_drop_const(X, "use", what = v)   # D-82
  X[, `:=`(rej = use & p < 0.05, rejn = use & p < 0.05, estimate = NA_real_, se = NA_real_)]
  ind[[v]] <- X
}
# ---- per design point: type I error (H0), power and usable share -------------------------------------------------------
ND <- BASE[, .(nd = sum(!design_degenerate)), by = .(scen, N, ntime, label)]
DP <- rbindlist(lapply(names(ind), function(s) ind[[s]][, .(n = sum(use), k = sum(rej), kn = sum(rejn)), by = .(scen, N, ntime, label)][, strategy := s]))
DP <- merge(DP, ND, by = c("scen", "N", "ntime", "label")); DP[, dgp := DGP_COP]
mk <- function(kcol, ref) { H0 <- flags(DP[label == "manual", .(dgp, scen, N, ntime, strategy, n, nd, kk = get(kcol))])
  H1 <- DP[label == "mmlb", .(dgp, scen, N, ntime, strategy, upow = 100 * get(kcol) / nd, pu_h1 = 100 * n / nd)]
  merge(H0, H1, by = c("dgp", "scen", "N", "ntime", "strategy"))[, ref := ref] }
R <- rbind(mk("k", "t"), mk("kn", "n")); R[, small := ntime <= 3 & N <= 25]
fwrite(R, file.path(TABDIR, "V4_S7e_dp_rates.csv"))
# ---- bias and coverage of the marginal estimators (both arms; true marginal value exact) ----------------------------------
MARG <- c("mmm_cc", "mmm_study", "gee_exch_unc", "gee_exch_md", "gee_exch_kc", "gee_exch_fg", "gee_ar1_unc", "gee_ar1_md", "gee_ar1_kc", "gee_ar1_fg", "gee_ind_fg", "shrink")
BC <- rbindlist(lapply(MARG, function(s) { X <- ind[[s]][use == TRUE]; X[, th := truth(scen, label)]
  X[, .(n = .N, bias = mean(estimate) - th[1], empse = sd(estimate), modse = sqrt(mean(se^2)),
        cov_n = 100 * mean(abs(estimate - th) <= qnorm(0.975) * se), cov_t = 100 * mean(abs(estimate - th) <= qt(0.975, N - np) * se)), by = .(scen, N, ntime, label)][, strategy := s] }))
BC[, `:=`(mcse_bias = empse / sqrt(n), bias_rel = bias / empse, mcse_cov_n = sqrt(cov_n * (100 - cov_n) / n))]
fwrite(BC, file.path(TABDIR, "V4_S7e_bias_cov.csv"))
# ---- usability reasons --------------------------------------------------------------------------------------------------
UR <- rbind(S[!(U0 %in% TRUE), .(why = substr(fifelse(is.na(status), "not usable (non-finite or |estimate| >= 10 or se >= 10)", status), 1, 70)), by = .(what = method)][, .N, by = .(what, why)],
            RECS[design_degenerate == FALSE & !is.finite(p_flip), .(what = "flip", why = flip_why)][, .N, by = .(what, why)])
fwrite(UR[order(what, -N)], file.path(TABDIR, "V4_S7e_usable_reasons.csv"))
# ---- rule of Figure 1 (estimators with estimate and interval; t at N - p) ---------------------------------------------------
EST <- c("mmm_cc", "glmer", "shrink", "gee_exch_unc", grep("^gee_(exch|ind|ar1)_(md|kc|fg|cr2|unc)$", names(ind), value = TRUE))
EST <- unique(EST)
RS <- fread(file.path(TABDIR, "V4_S7d_dp_rates.csv"))[ref == "t", .(dgp, scen, N, ntime, strategy, hi, lo, pu, rate, upow)]
RS[strategy == "study_mmm", strategy := "mmm_cc"]
RB <- fread(file.path(TABDIR, "V4_S7b_rule_dp.csv"))[ref == "t" & strategy %in% c("glmer|t", "shrink|t"), .(dgp, scen, N, ntime, strategy = sub("\\|t$", "", strategy), hi, lo, pu, rate, upow)]
OLD <- rbind(RS[strategy %in% EST], RB)
NEW <- R[ref == "t" & strategy %in% EST, .(dgp, scen, N, ntime, strategy, hi, lo, pu, rate, upow)]
rule <- function(D, npts, what) { x <- D[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), mean_pu = mean(pu), mean_rate = mean(rate, na.rm = TRUE), mean_upow = mean(upow)), by = strategy]
  x[, eligible := n_hi == 0L & n_dp == npts]; x[, chosen := eligible & mean_pu == max(mean_pu[eligible])]; x[, set := what]; setorder(x, -eligible, -mean_pu); x }
RU <- rbind(rule(OLD, 1896L, "study mechanisms (1,896 points)"), rule(NEW, 480L, "copula (480 points)"), rule(rbind(OLD, NEW), 2376L, "all four mechanisms (2,376 points)"))
print(RU[eligible == TRUE | strategy %in% c("mmm_cc", "gee_ar1_fg")], digits = 3); fwrite(RU, file.path(TABDIR, "V4_S7e_rule.csv"))
# ---- block bootstrap ---------------------------------------------------------------------------------------------------------
BI <- readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds"))
CNT <- lapply(1:5, function(m) { Mx <- BI[[m]]; Cm <- matrix(0, 2000L, 1000L); for (b in 1:2000) Cm[b, ] <- tabulate(Mx[b, ], 1000L); t(rbind(1, Cm)) })
DPK <- unique(R[, .(scen, N, ntime)])[order(scen, N, ntime)]; DPK[, key := paste(scen, N, ntime)]; DPK[, model := as.integer(sub("M", "", scen))]
NDI <- BASE[, .(scen, N, ntime, label, sim_id, nd = !design_degenerate)]; NDI[, key := paste(scen, N, ntime)]
# rates of an indicator table X (key, label, sim_id, use, rej): 480 x 2001 (column 1 = point estimate)
boot_mat <- function(X, arm = "manual", what = "rate", col = "rej") {
  X <- X[label == arm]; out <- matrix(NA_real_, nrow(DPK), 2001L)
  for (m in 1:5) {
    kk <- DPK[model == m]$key; x <- X[key %in% kk]
    U <- matrix(0, length(kk), 1000L); Kx <- U; ix <- cbind(match(x$key, kk), x$sim_id)
    U[ix] <- as.numeric(x$use); Kx[ix] <- as.numeric(x[[col]]); W <- CNT[[m]]
    if (what == "rate") out[match(kk, DPK$key), ] <- 100 * (Kx %*% W) / (U %*% W)
    else { y <- NDI[label == arm & key %in% kk]; Dn <- matrix(0, length(kk), 1000L); Dn[cbind(match(y$key, kk), y$sim_id)] <- as.numeric(y$nd)
           out[match(kk, DPK$key), ] <- 100 * (U %*% W) / (Dn %*% W) }
  }
  out
}
kx <- function(X) { X <- copy(X); X[, key := paste(scen, N, ntime)]; X }
summ <- function(M) { mm <- colMeans(M, na.rm = TRUE); c(est = mm[1], lo = unname(quantile(mm[-1], 0.025, na.rm = TRUE)), hi = unname(quantile(mm[-1], 0.975, na.rm = TRUE)), sd = sd(mm[-1], na.rm = TRUE), n_pts = sum(is.finite(M[, 1]))) }
# study mechanisms: flip and GEE AR(1) FG from the S7d bootstrap matrices (same indices), MMM cc and GLMM Wald from the study rows
BMS <- readRDS(file.path(REPO, "results_v4/S7d/agg/boot_mats.rds"))
old_mean <- function(M) { D <- copy(BMS$DPK); sel <- D$scen %in% sprintf("M%d", 1:5)
  A <- M[sel, , drop = FALSE]; D <- D[sel]; key <- paste(D$scen, D$N, D$ntime)
  out <- matrix(NA_real_, nrow(DPK), 2001L)
  for (k in seq_len(nrow(DPK))) { r <- which(key == DPK$key[k]); stopifnot(length(r) == 3L); out[k, ] <- colMeans(A[r, , drop = FALSE]) }
  out }
SR <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))[method %in% c("mmm", "glmer") & term %in% "time:trt" & scen %in% sprintf("M%d", 1:5),
        .(dgp, scen, N, ntime, label, method, sim_id, estimate, se, U0)]
SR <- merge(SR, NDI[, .(scen, N, ntime, label, sim_id, nd)], by = c("scen", "N", "ntime", "label", "sim_id"))
SR[, np := NP[scen]]; SR[, use := U0 & nd]; s9_drop_const(SR, "use", what = "study rows");   # D-82
s11_drop_cv(SR, "use", rows = SR$method == "mmm", what = "study rows")   # D-90
 SR[, rej := use & 2 * pt(-abs(estimate / se), N - np) < 0.05]; SR[, key := paste(scen, N, ntime)]
old_ind <- function(meth, what = "rate") { Ms <- lapply(unique(SR$dgp), function(g) boot_mat(SR[method == meth & dgp == g], what = what)); Reduce(`+`, Ms) / length(Ms) }
COP <- list(flip = boot_mat(kx(ind$flip)), ar1fg = boot_mat(kx(ind$gee_ar1_fg)), mmm = boot_mat(kx(ind$mmm_cc)), glmer = boot_mat(kx(ind$glmer)),
            mmm_pu = boot_mat(kx(ind$mmm_cc), what = "pu"))
OLDM <- list(flip = old_mean(BMS$FLIP), ar1fg = old_mean(BMS$ARFG), mmm = old_ind("mmm"), glmer = old_ind("glmer"), mmm_pu = old_ind("mmm", "pu"))
rm(SR); invisible(gc())
# bias / empSE and coverage of MMM cc under H1 (weighted sums with the bootstrap counts)
bias_boot <- function(X) { X <- kx(X)[label == "mmlb"]; out <- matrix(NA_real_, nrow(DPK), 2001L); cv <- out
  for (m in 1:5) { kk <- DPK[model == m]$key; x <- X[key %in% kk & use == TRUE]; th <- ifelse(m == 4L, 0.05, 0.2)
    E <- matrix(0, length(kk), 1000L); U <- E; C <- E; ix <- cbind(match(x$key, kk), x$sim_id)
    E[ix] <- x$estimate; U[ix] <- 1; C[ix] <- as.numeric(abs(x$estimate - th) <= qnorm(0.975) * x$se); W <- CNT[[m]]
    n <- U %*% W; s1 <- E %*% W; s2 <- (E^2) %*% W; mu <- s1 / n; sdv <- sqrt((s2 - n * mu^2) / (n - 1))
    out[match(kk, DPK$key), ] <- (mu - th) / sdv; cv[match(kk, DPK$key), ] <- 100 * (C %*% W) / n }
  list(rel = out, cov = cv) }
BB <- bias_boot(ind$mmm_cc)
BM <- rbindlist(list(
  data.table(q = "flip rate H0, copula", t(summ(COP$flip))), data.table(q = "flip rate H0, copula - study mechanisms", t(summ(COP$flip - OLDM$flip))),
  data.table(q = "GEE AR(1) FG t rate H0, copula", t(summ(COP$ar1fg))), data.table(q = "GEE AR(1) FG t rate H0, copula - study mechanisms", t(summ(COP$ar1fg - OLDM$ar1fg))),
  data.table(q = "MMM cc t rate H0, copula", t(summ(COP$mmm))), data.table(q = "MMM cc t rate H0, copula - study mechanisms", t(summ(COP$mmm - OLDM$mmm))),
  data.table(q = "GLMM Wald t rate H0, copula", t(summ(COP$glmer))), data.table(q = "GLMM Wald t rate H0, copula - study mechanisms", t(summ(COP$glmer - OLDM$glmer))),
  data.table(q = "MMM cc usable H0, copula", t(summ(COP$mmm_pu))), data.table(q = "MMM cc usable H0, copula - study mechanisms", t(summ(COP$mmm_pu - OLDM$mmm_pu))),
  data.table(q = "MMM cc bias/empSE H1, copula", t(summ(BB$rel))), data.table(q = "MMM cc coverage (normal) H1, copula", t(summ(BB$cov)))))
setnames(BM, c("q", "est", "lo", "hi", "sd", "n_pts")); print(BM, digits = 4); fwrite(BM, file.path(TABDIR, "V4_S7e_boot_means.csv"))
saveRDS(list(DPK = DPK, COP = COP, OLDM = OLDM, BB = BB), file.path(OUT, "agg/boot_mats.rds"))
# ---- predictions ----------------------------------------------------------------------------------------------------------
verdict <- function(lo, hi, a, b) fifelse(lo >= a & hi <= b, "confermata", fifelse(hi < a | lo > b, "smentita", "indeterminata"))
bm <- function(qq) { x <- BM[q == qq]; stopifnot(nrow(x) == 1L); x }
cnt <- function(s, col = "hi") R[ref == "t" & strategy == s, sum(get(col))]
row_bm <- function(code, qq, a, b) data.table(code = code, obs = bm(qq)$est, lo = bm(qq)$lo, hi = bm(qq)$hi, a = a, b = b)
row_ct <- function(code, v, a, b) data.table(code = code, obs = v, lo = v, hi = v, a = a, b = b)
ch <- RU[set == "all four mechanisms (2,376 points)" & chosen == TRUE]$strategy
P <- rbindlist(list(
  row_ct("PF1", cnt("flip"), 0, 0), row_bm("PF2", "flip rate H0, copula - study mechanisms", -1, 1),
  row_ct("PF3", cnt("gee_ar1_fg"), 0, 0), row_ct("PF4", cnt("mmm_cc"), 0, 0),
  row_bm("PF5", "GEE AR(1) FG t rate H0, copula - study mechanisms", -1, 1), row_bm("PF6", "MMM cc t rate H0, copula - study mechanisms", -1, 1),
  row_bm("PF7", "MMM cc bias/empSE H1, copula", -0.46, 0.46), row_bm("PF8", "MMM cc coverage (normal) H1, copula", 92.5, 97.5),
  row_bm("PF9", "GLMM Wald t rate H0, copula - study mechanisms", 1, 10),
  row_ct("PF10", as.numeric(identical(ch, "gee_ar1_fg")), 1, 1),
  row_bm("PF11", "MMM cc usable H0, copula - study mechanisms", -2.5, 2.5)))
P[, verdict := verdict(lo, hi, a, b)]
P[code == "PF10", extra := paste("chosen:", paste(ch, collapse = ", "))]
P[code == "PF7", extra := sprintf("H1 points with bias/empSE +/- 2 MCSE wholly outside +/- 0.46: %d", BC[strategy == "mmm_cc" & label == "mmlb", sum(abs(bias_rel) - 2 / sqrt(n) > 0.46, na.rm = TRUE)])]
print(P, digits = 4); fwrite(P, file.path(TABDIR, "V4_S7e_predictions.csv"))
tm <- colMeans(SECS, na.rm = TRUE); cat("mean seconds per replicate by part:", paste(names(tm), sprintf("%.3f", tm), collapse = " "), "\n")
cat(sprintf("done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
