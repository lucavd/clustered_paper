# =============================================================================
# s7d_agg.R -- S7d aggregation: per design point x arm x strategy rates with MCSE and Bradley flags, the rule of
# Figure 1 reapplied with the new strategies, block-bootstrap aggregates (B = 2000, indices of S4 and S6b), and the
# pre-registered predictions PE1-PE12 with their verdicts. Gate GE1 (completeness) is computed here too.
# Inputs: results_v4/S7d/rows/*.rds; study rows results_v4/variant_cc/S4/rows_v4.rds, results_v4/S6b/ns/rows_ns.rds;
#         tables/V4_S7b_rule_dp.csv (study strategies, t at N - p).
# Outputs: tables/V4_S7d_dp_rates.csv, V4_S7d_rule.csv, V4_S7d_predictions.csv, V4_S7d_GE1.csv, V4_S7d_boot_means.csv,
#          V4_S7d_usable_reasons.csv; results_v4/S7d/agg/ (per-replicate indicator matrices).
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R"))       # REPO, TABDIR, NP, WALD7, flags()
suppressMessages(library(data.table)); setDTthreads(8L)
t0 <- Sys.time(); OUT <- file.path(REPO, "results_v4/S7d"); dir.create(file.path(OUT, "agg"), showWarnings = FALSE)
source(path.expand("~/MMM/R/s7d/s7d_common.R"))
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
source(path.expand("~/MMM/R/s11/s11_cv.R"))   # D-90 (S11, D-91): no result for the MMM with the check when its refit does not converge
TK <- s7d_tasks()
# ---- read ---------------------------------------------------------------------------------------------------
have <- file.exists(TK$file)
L <- parallel::mclapply(TK$file[have], readRDS, mc.cores = 16L)
RECS <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); ROWS <- rbindlist(lapply(L, `[[`, "rows"), fill = TRUE); rm(L); invisible(gc())
cat(sprintf("read %d chunks: %d recs, %d rows (%.1f min)\n", sum(have), nrow(RECS), nrow(ROWS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
# ---- GE1 completeness ---------------------------------------------------------------------------------------
nrow_per <- ROWS[, .(cnt = .N), by = K6]
gen_bad_nondeg <- RECS[generated == FALSE & design_degenerate == FALSE, .N]
ge1 <- data.table(chunks = sum(have), chunks_expected = nrow(TK), recs = nrow(RECS), recs_unique = uniqueN(RECS[, ..K6]),
                  generated = RECS[generated == TRUE, .N], gen_errors = RECS[generated == FALSE, .N], gen_errors_nondegenerate = gen_bad_nondeg,
                  rows_per_generated_ok = all(nrow_per$cnt == 3L) && nrow(nrow_per) == RECS[generated == TRUE, .N])
ge1[, pass_completeness := chunks == chunks_expected & recs == 3792000L & recs_unique == recs & gen_errors_nondegenerate == 0L & rows_per_generated_ok]
print(ge1); fwrite(ge1, file.path(TABDIR, "V4_S7d_GE1.csv"))
# ---- per-replicate indicators of the new strategies ---------------------------------------------------------
ROWS <- merge(ROWS, RECS[, c(K6, "design_degenerate"), with = FALSE], by = K6)
ROWS[, np := NP[scen]]
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
ind_list <- list()
for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) {
  if (ty == "cr2" && cn == "ar1") next
  X <- ROWS[cor == cn, c(K6, "est", "np", "design_degenerate", "gee_ok", SEC[[ty]]), with = FALSE]; setnames(X, SEC[[ty]], "se")
  X[, use := !design_degenerate & gee_ok & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(X, "use", what = paste0("gee_", cn, "_", ty))   # D-82
  X[, z := est / se]
  X[, `:=`(rej = use & 2 * pt(-abs(z), N - np) < 0.05, rej2 = use & 2 * pt(-abs(z), N - 2L) < 0.05, rejn = use & 2 * pnorm(-abs(z)) < 0.05)]
  ind_list[[paste0("gee_", cn, "_", ty)]] <- X[, c(K6, "use", "rej", "rej2", "rejn"), with = FALSE]
}
for (v in c("flip", "flip_strict")) {
  pc <- if (v == "flip") "p_flip" else "p_flip_strict"
  X <- RECS[, c(K6, "design_degenerate", pc), with = FALSE]; setnames(X, pc, "p")
  X[, use := !design_degenerate & is.finite(p)]; s9_drop_const(X, "use", what = v); X[, rej := use & p < 0.05]; X[, `:=`(rej2 = rej, rejn = rej)]   # D-82
  ind_list[[v]] <- X[, c(K6, "use", "rej", "rej2", "rejn"), with = FALSE]
}
ind_list[["nondeg"]] <- RECS[, c(K6, "design_degenerate"), with = FALSE][, `:=`(use = !design_degenerate, rej = FALSE, rej2 = FALSE, rejn = FALSE)][, c(K6, "use", "rej", "rej2", "rejn"), with = FALSE]
# study strategies needed per replicate (bootstrap of paired comparisons): saved rows, U0 filter, same denominators
KEEP <- c("dgp", "scen", "N", "ntime", "label", "method", "sim_id", "term", "estimate", "se", "U0")
SA <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))[method %in% c("gee", "gee_md", "gee_kc", "gee_fg", "mmm") & term %in% "time:trt", ..KEEP]
SB <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds"))[method %in% c("gee", "gee_md", "gee_kc", "gee_fg", "mmm") & term %in% "time:trt", ..KEEP]
S <- rbind(SA, SB); rm(SA, SB); invisible(gc())
S <- merge(S, RECS[, c(K6, "design_degenerate"), with = FALSE], by = K6)
S[, np := NP[scen]]; S[, use := U0 & !design_degenerate]; s9_drop_const(S, "use", what = "study rows"); S[, z := estimate / se]   # D-82
s11_drop_cv(S, "use", rows = S$method == "mmm", what = "study rows")   # D-90
S[, `:=`(rej = use & 2 * pt(-abs(z), N - np) < 0.05, rej2 = use & 2 * pt(-abs(z), N - 2L) < 0.05, rejn = use & 2 * pnorm(-abs(z)) < 0.05)]
for (m in unique(S$method)) ind_list[[paste0("study_", m)]] <- S[method == m, c(K6, "use", "rej", "rej2", "rejn"), with = FALSE]
rm(S); invisible(gc())
# ---- per design point -----------------------------------------------------------------------------------------
ND <- RECS[, .(nd = sum(!design_degenerate)), by = .(dgp, scen, N, ntime, label)]
DP <- rbindlist(lapply(setdiff(names(ind_list), "nondeg"), function(s) ind_list[[s]][, .(n = sum(use), k = sum(rej), k2 = sum(rej2), kn = sum(rejn)), by = .(dgp, scen, N, ntime, label)][, strategy := s]))
DP <- merge(DP, ND, by = c("dgp", "scen", "N", "ntime", "label"))
DP[, arm := fifelse(label == "mmlb", "H1", "H0")]
mk <- function(kcol, ref) { H0 <- flags(DP[arm == "H0", .(dgp, scen, N, ntime, strategy, n, nd, kk = get(kcol))])
  H1 <- DP[arm == "H1", .(dgp, scen, N, ntime, strategy, upow = 100 * get(kcol) / nd, pu_h1 = 100 * n / nd)]
  merge(H0, H1, by = c("dgp", "scen", "N", "ntime", "strategy"))[, ref := ref] }
R <- rbind(mk("k", "t"), mk("k2", "t2"), mk("kn", "n"))
R[, small := ntime <= 3 & N <= 25]
fwrite(R, file.path(TABDIR, "V4_S7d_dp_rates.csv"))
# ---- usability reasons (new strategies) ---------------------------------------------------------------------------
UR <- rbind(ROWS[!design_degenerate & !gee_ok, .(why = substr(gee_msg, 1, 60), what = paste0("geeglm_", cor))][, .N, by = .(what, why)],
            rbindlist(lapply(c("md", "kc", "fg", "cr2"), function(ty) ROWS[!design_degenerate & gee_ok & !is.na(get(paste0("msg_", ty))),
              .(what = paste0(ty, "_", cor), why = substr(get(paste0("msg_", ty)), 1, 60))][, .N, by = .(what, why)])),
            RECS[!design_degenerate & !is.finite(p_flip), .(what = "flip", why = flip_why)][, .N, by = .(what, why)])
fwrite(UR[order(what, -N)], file.path(TABDIR, "V4_S7d_usable_reasons.csv"))
# ---- rule of Figure 1 -----------------------------------------------------------------------------------------------
ST <- fread(file.path(TABDIR, "V4_S7b_rule_dp.csv"))[ref == "t", .(dgp, scen, N, ntime, strategy = paste0("paper_", sub("\\|t$", "", strategy)), n, nd, rate, mc, assess, hi, lo, pu, upow)]
ST[, small := ntime <= 3 & N <= 25]
RU <- rbind(ST[, .(src = "study (saved rows)", n_dp = .N, mean_rate = mean(rate, na.rm = TRUE), n_hi = sum(hi), n_lo = sum(lo), mean_pu = mean(pu), mean_upow = mean(upow),
                   pu_small = mean(pu[small]), rate_small = mean(rate[small], na.rm = TRUE), hi_small = sum(hi[small])), by = strategy],
            R[ref == "t", .(src = "S7d run", n_dp = .N, mean_rate = mean(rate, na.rm = TRUE), n_hi = sum(hi), n_lo = sum(lo), mean_pu = mean(pu), mean_upow = mean(upow),
                   pu_small = mean(pu[small]), rate_small = mean(rate[small], na.rm = TRUE), hi_small = sum(hi[small])), by = strategy])
RU[, wald := !startsWith(strategy, "flip")]
RU[, eligible := n_hi == 0L & n_dp == 1896L]
RU[, chosen_wald := wald & eligible & mean_pu == max(mean_pu[wald & eligible])]
RU[, chosen_all := eligible & mean_pu == max(mean_pu[eligible])]
setorder(RU, -eligible, -mean_pu); print(RU, digits = 3)
fwrite(RU, file.path(TABDIR, "V4_S7d_rule.csv"))
# ---- block bootstrap ------------------------------------------------------------------------------------------------
BI <- c(readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds")), readRDS(file.path(REPO, "results_v4/S6b/ns/boot_idx_67.rds"))); names(BI) <- as.character(1:7)
CNT <- lapply(BI, function(Mx) { Cm <- matrix(0, 2000L, 1000L); for (b in 1:2000) Cm[b, ] <- tabulate(Mx[b, ], 1000L); t(rbind(1, Cm)) })   # 1000 x 2001
model_of <- function(scen) fifelse(scen == "M1c", 1L, as.integer(sub("^M", "", sub("c$", "", scen))))
DPK <- unique(R[, .(dgp, scen, N, ntime)])[order(dgp, scen, N, ntime)]; DPK[, key := paste(dgp, scen, N, ntime)]; DPK[, model := model_of(scen)]
ND1 <- ind_list[["nondeg"]][label == "manual"]; ND1[, key := paste(dgp, scen, N, ntime)]
boot_rates <- function(s, arm = "H0", what = "rate") {   # 1896 x 2001 matrix (column 1 = point estimate)
  X <- ind_list[[s]][label == (if (arm == "H0") "manual" else "mmlb")]; X[, key := paste(dgp, scen, N, ntime)]
  out <- matrix(NA_real_, nrow(DPK), 2001L)
  for (m in unique(DPK$model)) {
    kk <- DPK[model == m]$key; x <- X[key %in% kk]
    U <- matrix(0, length(kk), 1000L); Kx <- U; ix <- cbind(match(x$key, kk), x$sim_id)
    U[ix] <- as.numeric(x$use); Kx[ix] <- as.numeric(x$rej)
    W <- CNT[[as.character(m)]]
    if (what == "rate") out[match(kk, DPK$key), ] <- 100 * (Kx %*% W) / (U %*% W)
    else { y <- ND1[key %in% kk]; Dn <- matrix(0, length(kk), 1000L); Dn[cbind(match(y$key, kk), y$sim_id)] <- as.numeric(y$use)
           out[match(kk, DPK$key), ] <- 100 * (U %*% W) / (Dn %*% W) }
  }
  out
}
# what = "pu": usable share over the replicates with an estimable design (denominator of the study, D-11)
summ <- function(M, sel = rep(TRUE, nrow(M))) { mm <- colMeans(M[sel, , drop = FALSE], na.rm = TRUE)
  c(est = mm[1], lo = unname(quantile(mm[-1], 0.025)), hi = unname(quantile(mm[-1], 0.975)), sd = sd(mm[-1])) }
FLIP <- boot_rates("flip"); FCFG <- boot_rates("gee_exch_fg"); STFG <- boot_rates("study_gee_fg")
INDFG <- boot_rates("gee_ind_fg"); ARFG <- boot_rates("gee_ar1_fg"); FCKC <- boot_rates("gee_exch_kc"); STKC <- boot_rates("study_gee_kc")
PU_FCFG <- boot_rates("gee_exch_fg", what = "pu"); PU_STFG <- boot_rates("study_gee_fg", what = "pu")
small <- DPK$N <= 25 & DPK$ntime <= 3
BM <- rbindlist(list(
  data.table(q = "flip rate H0", t(summ(FLIP))), data.table(q = "flip_strict rate H0", t(summ(boot_rates("flip_strict")))),
  data.table(q = "fewclust FG exch - study FG rate", t(summ(FCFG - STFG))),
  data.table(q = "FG ind - FG exch rate (fewclust)", t(summ(INDFG - FCFG))), data.table(q = "FG ar1 - FG exch rate (fewclust)", t(summ(ARFG - FCFG))),
  data.table(q = "fewclust KC exch rate", t(summ(FCKC))), data.table(q = "study KC (geesmv) rate", t(summ(STKC))), data.table(q = "fewclust KC - study KC rate", t(summ(FCKC - STKC))),
  data.table(q = "pu fewclust FG exch - study FG", t(summ(PU_FCFG - PU_STFG))),
  data.table(q = "fewclust FG exch rate", t(summ(FCFG))), data.table(q = "study FG rate", t(summ(STFG)))))
setnames(BM, c("q", "est", "lo", "hi", "sd")); print(BM, digits = 4); fwrite(BM, file.path(TABDIR, "V4_S7d_boot_means.csv"))
saveRDS(list(DPK = DPK, FLIP = FLIP, FCFG = FCFG, STFG = STFG, INDFG = INDFG, ARFG = ARFG, FCKC = FCKC, STKC = STKC), file.path(OUT, "agg/boot_mats.rds"))
# ---- predictions ------------------------------------------------------------------------------------------------------
Rt <- R[ref == "t"]; W <- dcast(Rt, dgp + scen + N + ntime + small ~ strategy, value.var = c("pu", "hi", "lo", "rate"))
STw <- dcast(ST, dgp + scen + N + ntime ~ strategy, value.var = c("pu", "hi", "rate"))
W <- merge(W, STw, by = c("dgp", "scen", "N", "ntime"))
verdict <- function(lo, hi, a, b) fifelse(lo >= a & hi <= b, "confermata", fifelse(hi < a | lo > b, "smentita", "indeterminata"))
bm <- function(qq) { x <- BM[q == qq]; stopifnot(nrow(x) == 1L); x }
new_wald <- RU[src == "S7d run" & wald & !startsWith(strategy, "study_") & strategy != "gee_exch_unc"]
P <- rbindlist(list(
  data.table(code = "PE1", obs = bm("flip rate H0")$est, lo = bm("flip rate H0")$lo, hi = bm("flip rate H0")$hi, a = 2.5, b = 7.5),
  data.table(code = "PE2", obs = W[, sum(hi_flip)], lo = W[, sum(hi_flip)], hi = W[, sum(hi_flip)], a = 0, b = 0),
  data.table(code = "PE3", obs = W[, sum(lo_flip)], lo = W[, sum(lo_flip)], hi = W[, sum(lo_flip)], a = 0, b = 60),
  data.table(code = "PE4", obs = W[, min(pu_flip - pu_paper_gee_fg)], lo = W[, min(pu_flip - pu_paper_gee_fg)], hi = Inf, a = -2.5, b = Inf),
  data.table(code = "PE5", obs = W[small == TRUE, mean(pu_flip)], lo = W[small == TRUE, mean(pu_flip)], hi = W[small == TRUE, mean(pu_flip)], a = 90, b = 100),
  data.table(code = "PE6", obs = bm("fewclust FG exch - study FG rate")$est, lo = bm("fewclust FG exch - study FG rate")$lo, hi = bm("fewclust FG exch - study FG rate")$hi, a = -1, b = 1),
  data.table(code = "PE7", obs = W[, sum(hi_gee_exch_fg != hi_paper_gee_fg)], lo = W[, sum(hi_gee_exch_fg != hi_paper_gee_fg)], hi = W[, sum(hi_gee_exch_fg != hi_paper_gee_fg)], a = 0, b = 2),
  data.table(code = "PE8", obs = bm("pu fewclust FG exch - study FG")$est, lo = bm("pu fewclust FG exch - study FG")$lo, hi = bm("pu fewclust FG exch - study FG")$hi, a = 0, b = 2),
  data.table(code = "PE9", obs = bm("FG ind - FG exch rate (fewclust)")$est, lo = bm("FG ind - FG exch rate (fewclust)")$lo, hi = bm("FG ind - FG exch rate (fewclust)")$hi, a = -1, b = 1),
  data.table(code = "PE10", obs = bm("FG ar1 - FG exch rate (fewclust)")$est, lo = bm("FG ar1 - FG exch rate (fewclust)")$lo, hi = bm("FG ar1 - FG exch rate (fewclust)")$hi, a = -1, b = 1),
  data.table(code = "PE11", obs = bm("fewclust KC exch rate")$est, lo = bm("fewclust KC exch rate")$lo, hi = bm("fewclust KC exch rate")$hi, a = 4.8, b = 6.8),
  data.table(code = "PE12", obs = new_wald[, sum(eligible)], lo = new_wald[, sum(eligible)], hi = new_wald[, sum(eligible)], a = 0, b = 0)))
P[, verdict := verdict(lo, hi, a, b)]
P[code == "PE8", extra := sprintf("min per point %.2f", W[, min(pu_gee_exch_fg - pu_paper_gee_fg)])]
P[code == "PE5", extra := sprintf("rate small %.2f, hi small %d", W[small == TRUE, mean(rate_flip, na.rm = TRUE)], W[small == TRUE, sum(hi_flip)])]
print(P, digits = 4); fwrite(P, file.path(TABDIR, "V4_S7d_predictions.csv"))
cat(sprintf("done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
