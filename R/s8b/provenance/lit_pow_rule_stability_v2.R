#!/usr/bin/env Rscript
# =============================================================================
# lit_pow_rule_stability_v2.R -- S8, track Letteratura, giro 2: PoW-1 with the project's block-bootstrap indices (read-only on ~/MMM).
# v2 changes: (1) resampling by replicate blocks (model, i) with the stored indices of the project
#   (results_v4/S4/boot_idx.rds for models 1-5, results_v4/S6b/ns/boot_idx_67.rds for 6-7; M1c -> model 1;
#   Copula-AR1 scenario Mk -> model k), as in R/s7d/s7d_agg.R: the same index draw b is applied to every point of a
#   model, across mechanisms, N and T (common random numbers); (2) plug-in with empirical-Bayes shrunken rates
#   (within strategy x mechanism) as an optimistic counterpart of the raw plug-in; (3) outputs suffixed _v2.
# PoW-1: stability of the decision rule of Section 2.6 (Table 4) under resampling of replicate indices.
#   Rows: results_v4/S7d/rows (1,896 points: main grid, M1c, M6-M7) and results_v4/S7e/rows (480 Copula-AR1 points),
#   H0 arm, GEE with exchangeable / independence / AR(1) working correlation and the fewclust corrections (MD, KC, FG)
#   plus the Ford-Westgate average of MD and KC (average of the two variances), Wald t with N - np df (np as in
#   R/s7b/s7b_common.R). MMM with the check from results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds (its GEE AR(1) FG
#   rows are a cross-check). Usable = !design_degenerate & gee_ok & |est| < 10 & 0 < se < 10 (R/s7d/s7d_agg.R).
#   Gate: per-point n, k, nd equal tables/V5_S7f_cp_S7d_dp_rates.csv and _S7e_ (ref t); flag counts equal
#   tables/V5_S7f_table4.csv.
#   Bootstrap: B = 2000, seed 20261015, replicate indices 1..1000 resampled jointly across all cells (Section 2.5).
#   Per resample: per point n*, k*, nd*; "clearly above 7.5%" = n* >= 20 & Clopper-Pearson lower limit > 7.5%;
#   eligible = no flag; the rule picks the eligible strategy with the largest mean usable share.
#   Also the plug-in binomial probability of zero flags (points independent, true rate = observed rate).
# PoW-3: why GEE AR(1) FG returns no usable result with T <= 3 and N <= 25 (both arms).
# Run: cd <workdir> && OMP_NUM_THREADS=1 Rscript --vanilla lit_pow_rule_stability.R
# =============================================================================
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(1L)
t0 <- Sys.time()
REPO <- "/home/user/MMM"
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
cp_lo <- function(k, n) ifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
cp_up <- function(k, n) ifelse(k >= n, 1, qbeta(0.975, k + 1, pmax(n - k, 1e-12)))
files <- c(list.files(file.path(REPO, "results_v4/S7d/rows"), full.names = TRUE),
           list.files(file.path(REPO, "results_v4/S7e/rows"), full.names = TRUE))
L <- parallel::mclapply(files, function(f) { x <- readRDS(f)
  g <- as.data.table(if (!is.null(x$gee)) x$gee else x$rows); r <- as.data.table(x$recs)
  list(recs = r[, .(sim_id, design_degenerate, flip_why, sum_y, n_obs, dgp, scen, N, ntime, label)],
       gee = g[, .(cor, gee_ok, est, se_md, se_kc, se_fg, msg_fg, sim_id, dgp, scen, N, ntime, label)]) }, mc.cores = 2L)
RECS <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); G <- rbindlist(lapply(L, `[[`, "gee"), fill = TRUE); rm(L); invisible(gc())
cat(sprintf("read %d files: %d recs, %d gee rows (%.1f min)\n", length(files), nrow(RECS), nrow(G), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
G <- merge(G, RECS[, c(K6, "design_degenerate", "flip_why"), with = FALSE], by = K6)
G[, np := NP[scen]]
G[, se_av := sqrt((se_md^2 + se_kc^2) / 2)]
SEC <- c(fg = "se_fg", md = "se_md", kc = "se_kc", avg = "se_av")
IND <- list()
for (cn in c("ar1", "ind", "exch")) for (ty in names(SEC)) {
  X <- G[cor == cn, c(K6, "est", "np", "design_degenerate", "gee_ok", SEC[[ty]]), with = FALSE]; setnames(X, SEC[[ty]], "se")
  X[, use := !design_degenerate & gee_ok & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(X, "use", what = paste0("gee_", cn, "_", ty))   # D-82
  X[, rej := use & 2 * pt(-abs(est / se), N - np) < 0.05]
  IND[[paste0("gee_", cn, "_", ty)]] <- X[, c(K6, "use", "rej"), with = FALSE]
}
PR <- as.data.table(readRDS(file.path(REPO, "results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds")))
IND[["mmm_cc"]] <- PR[strategy == "mmm_cc", c(K6, "use", "rej"), with = FALSE]
IND[["chk_gee_ar1_fg"]] <- PR[strategy == "gee_ar1_fg", c(K6, "use", "rej"), with = FALSE]
ND <- RECS[, .(dgp, scen, N, ntime, label, sim_id, nondeg = !design_degenerate)]
pts <- unique(ND[label == "manual", .(dgp, scen, N, ntime)])[order(dgp, scen, N, ntime)]
pts[, pid := .I]; npt <- nrow(pts); cat("points:", npt, "\n")
A <- rbind(fread(file.path(REPO, "tables/V5_S7f_cp_S7d_dp_rates.csv")), fread(file.path(REPO, "tables/V5_S7f_cp_S7e_dp_rates.csv")))[ref == "t"]
A[strategy == "study_mmm", strategy := "mmm_cc"]
mat <- function(D, col) { M <- matrix(0, npt, 1000); D <- merge(D, pts, by = c("dgp", "scen", "N", "ntime"))
  M[cbind(D$pid, D$sim_id)] <- as.numeric(D[[col]]); M }
NDm <- mat(ND[label == "manual"], "nondeg")
USE <- list(); REJ <- list(); GT <- list()
for (s in names(IND)) {
  D <- IND[[s]][label == "manual"]
  USE[[s]] <- mat(D, "use"); REJ[[s]] <- mat(D, "rej")
  ref_s <- if (s == "chk_gee_ar1_fg") "gee_ar1_fg" else s
  R <- A[strategy == ref_s, .(dgp, scen, N, ntime, n_ref = n, k_ref = kk, nd_ref = nd)]
  P <- merge(pts, R, by = c("dgp", "scen", "N", "ntime"), all.x = TRUE)[order(pid)]
  GT[[s]] <- data.table(strategy = s, points_in_ref = sum(!is.na(P$n_ref)), n_diff = sum(rowSums(USE[[s]]) != P$n_ref, na.rm = TRUE),
                        k_diff = sum(rowSums(REJ[[s]]) != P$k_ref, na.rm = TRUE), nd_diff = sum(rowSums(NDm) != P$nd_ref, na.rm = TRUE))
}
GT <- rbindlist(GT)
flag_hi <- function(k, n) n >= 20 & 100 * cp_lo(k, n) > 7.5
flag_lo <- function(k, n) n >= 20 & 100 * cp_up(k, n) < 2.5
kstar <- function(n) sapply(n, function(m) { if (m < 20) return(m + 1); kk <- 0:m; w <- which(100 * cp_lo(kk, m) > 7.5); if (length(w)) kk[min(w)] else m + 1 })
OBS <- rbindlist(lapply(names(IND), function(s) { n <- rowSums(USE[[s]]); k <- rowSums(REJ[[s]]); ks <- kstar(n)
  pf <- ifelse(n >= 20, pbinom(ks - 1, n, k / pmax(n, 1), lower.tail = FALSE), 0)
  data.table(strategy = s, hi_obs = sum(flag_hi(k, n)), lo_obs = sum(flag_lo(k, n)), pu_obs = mean(100 * n / rowSums(NDm)),
             max_rate = max((100 * k / n)[n >= 20]), plugin_exp_flags = sum(pf), plugin_p_pass = exp(sum(log1p(-pmin(pf, 1 - 1e-15))))) }))

EB <- rbindlist(lapply(names(IND), function(s) { n <- rowSums(USE[[s]]); k <- rowSums(REJ[[s]]); ph <- ifelse(n > 0, k / n, NA_real_)
  ps <- ph
  for (g in unique(pts$dgp)) { w <- which(pts$dgp == g & n >= 20); m <- sum(k[w]) / sum(n[w]); v <- m * (1 - m) / n[w]
    tau2 <- max(0, var(ph[w]) - mean(v)); ps[w] <- m + tau2 / (tau2 + v) * (ph[w] - m) }
  ks <- kstar(n); pf <- ifelse(n >= 20, pbinom(ks - 1, n, pmin(pmax(ps, 0), 1), lower.tail = FALSE), 0)
  data.table(strategy = s, eb_exp_flags = sum(pf, na.rm = TRUE), eb_p_pass = exp(sum(log1p(-pmin(pf, 1 - 1e-15)), na.rm = TRUE))) }))
OBS <- merge(OBS, EB, by = "strategy")
T4 <- fread(file.path(REPO, "tables/V5_S7f_table4.csv"))
GT <- merge(GT, OBS[, .(strategy, hi_obs, lo_obs, pu_obs)], by = "strategy", all.x = TRUE)
GT <- merge(GT, T4[, .(strategy, hi_all_t4 = hi_all, lo_all_t4 = lo_all, pu_t4 = pu)], by = "strategy", all.x = TRUE)
fwrite(GT, "lit_pow1_gate_v2.csv"); print(GT)
BI <- c(readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds")), readRDS(file.path(REPO, "results_v4/S6b/ns/boot_idx_67.rds"))); names(BI) <- as.character(1:7)
B <- 2000L
CNT <- lapply(BI, function(Mx) { Cm <- matrix(0, 1000L, B); for (b in 1:B) Cm[, b] <- tabulate(Mx[b, ], 1000L); Cm })
model_of <- function(scen) ifelse(scen == "M1c", 1L, as.integer(sub("^M", "", scen)))
pts[, model := model_of(scen)]
bprod <- function(M) { out <- matrix(NA_real_, npt, B); for (m in unique(pts$model)) { r <- pts$pid[pts$model == m]; out[r, ] <- M[r, , drop = FALSE] %*% CNT[[as.character(m)]] }; out }
NDb <- bprod(NDm)
strat_rule <- setdiff(names(IND), "chk_gee_ar1_fg")
HI <- matrix(NA_real_, B, length(strat_rule), dimnames = list(NULL, strat_rule)); PU <- HI; LO <- HI
for (s in strat_rule) { nb <- bprod(USE[[s]]); kb <- bprod(REJ[[s]])
  HI[, s] <- colSums(flag_hi(kb, nb)); LO[, s] <- colSums(flag_lo(kb, nb)); PU[, s] <- colMeans(100 * nb / NDb) }
PASS <- data.table(strategy = strat_rule, p_pass_boot = colMeans(HI == 0), mean_flags_boot = colMeans(HI),
                   q90_flags_boot = apply(HI, 2, quantile, 0.9), mean_lo_boot = colMeans(LO), mean_pu_boot = colMeans(PU))
PASS <- merge(PASS, OBS, by = "strategy"); setorder(PASS, -p_pass_boot)
fwrite(PASS, "lit_pow1_pass_v2.csv"); print(PASS, digits = 3)
ch <- sapply(seq_len(B), function(b) { e <- HI[b, ] == 0; if (!any(e)) "none" else strat_rule[e][which.max(PU[b, e])] })
CH <- data.table(choice = ch)[, .N, by = choice][, share := N / B][order(-N)]
fwrite(CH, "lit_pow1_choice_v2.csv"); print(CH)
JO <- data.table(event = c("FG AR(1) and MMM both eligible", "FG independence eligible", "FG exchangeable eligible",
                           "FG AR(1) eligible and FG independence not", "FG independence eligible and FG AR(1) not",
                           "pu FG independence > pu FG AR(1)"),
                 share = c(mean(HI[, "gee_ar1_fg"] == 0 & HI[, "mmm_cc"] == 0), mean(HI[, "gee_ind_fg"] == 0), mean(HI[, "gee_exch_fg"] == 0),
                           mean(HI[, "gee_ar1_fg"] == 0 & HI[, "gee_ind_fg"] > 0), mean(HI[, "gee_ind_fg"] == 0 & HI[, "gee_ar1_fg"] > 0),
                           mean(PU[, "gee_ind_fg"] > PU[, "gee_ar1_fg"])))
fwrite(JO, "lit_pow1_joint_v2.csv"); print(JO)
FP <- rbindlist(lapply(c("gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "mmm_cc", "gee_ar1_avg", "gee_exch_avg", "gee_ind_avg"), function(s) {
  n <- rowSums(USE[[s]]); k <- rowSums(REJ[[s]]); lo <- 100 * cp_lo(k, n); up <- 100 * cp_up(k, n)
  cbind(pts, strategy = s, n = n, k = k, rate = 100 * k / n, cp_lo = lo, cp_up = up)[n >= 20 & (lo > 6.5)] }))
fwrite(FP[order(strategy, -cp_lo)], "lit_pow1_flagged_points_v2.csv")
X <- G[cor == "ar1" & ntime <= 3 & N <= 25 & !design_degenerate]
X[, reason := fifelse(!gee_ok, "geeglm error",
              fifelse(!is.finite(est) | abs(est) >= 10, "estimate not finite or |est| >= 10",
              fifelse(!is.finite(se_fg), fifelse(grepl("singular", msg_fg), "FG SE: singular information (possible separation)", "FG SE: other failure"),
              fifelse(se_fg <= 0 | se_fg >= 10, "FG SE outside (0, 10)", "usable"))))]
X[, flip_nc := !is.na(flip_why)]
R3 <- X[, .N, by = .(label, reason, flip_not_computable = flip_nc)][, share := N / sum(N), by = label][order(label, -N)]
fwrite(R3, "lit_pow3_reasons_v2.csv"); print(R3)
R3b <- X[, .(n = .N, usable = mean(reason == "usable"), flip_nc = mean(flip_nc), est_bad = mean(reason == "estimate not finite or |est| >= 10"),
             se_sing = mean(grepl("singular", reason)), geeglm_err = mean(reason == "geeglm error")), by = .(label, ntime)]
fwrite(R3b, "lit_pow3_by_T_v2.csv"); print(R3b)
R3c <- X[, .(design_points = uniqueN(paste(dgp, scen, N, ntime)), datasets = .N, usable = sum(reason == "usable"),
             est_bad = sum(reason == "estimate not finite or |est| >= 10"), se_singular = sum(grepl("singular", reason)),
             se_other = sum(reason %in% c("FG SE: other failure", "FG SE outside (0, 10)")), geeglm_err = sum(reason == "geeglm error")), by = label]
R3c[, `:=`(share_sep_signature = (est_bad + se_singular) / datasets, share_usable = usable / datasets)]
fwrite(R3c, "lit_pow3_totals_v2.csv"); print(R3c)
writeLines(capture.output(sessionInfo()), "lit_pow_sessionInfo_v2.txt")
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
