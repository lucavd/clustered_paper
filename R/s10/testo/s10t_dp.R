# s10t_dp.R -- S10 "Testo": per-replicate usability and rejections written from the definitions of Sections 2.4-2.5
# of the manuscript, aggregated per design point x arm x strategy. Written from scratch; no project aggregation
# function is sourced. Inputs: results_v4/S10/testo/x_*.qs2 (copies made by s10t_extract*.R),
# tables/V4_S2_truth_glmmnative.csv (true marginal values, GLMM-native, for coverage only).
# Outputs: results_v4/S10/testo/dp_agg.csv.gz, cov_agg.csv.gz, cc_agg.csv.gz, ind_boot.qs2, checks.txt
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(3L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
CHK <- file(file.path(OUT, "checks.txt"), "w"); ck <- function(...) { s <- sprintf(...); cat(s, "\n"); writeLines(s, CHK) }
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)   # mean-model coefficients (Section 2.4: N - p)
# ---- base: every dataset of the 2,376 design points x 2 arms x 1,000 replicates ---------------------------------
R <- qs_read(file.path(OUT, "x_recs.qs2"), nthreads = 2L)
ck("recs %d unique keys %d", nrow(R), uniqueN(R[, ..K6]))
R[, est_design := generated & !design_degenerate]
R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
ck("generated FALSE: %d (of which non-degenerate %d); degenerate %d; constant outcome among estimable %d; constant among all generated %d",
   R[generated == FALSE, .N], R[generated == FALSE & !design_degenerate, .N], R[design_degenerate == TRUE, .N], R[const == TRUE, .N],
   R[generated == TRUE & (sum_y == 0 | sum_y == n_obs), .N])
CK <- as.data.table(readRDS("results_v4/S9/const_keys.rds")); CK[, inkeys := TRUE]
z <- merge(R[const == TRUE, ..K6][, mine := TRUE], CK, by = K6, all = TRUE)
ck("const: mine %d, project keys %d, both %d", z[mine %in% TRUE, .N], z[inkeys %in% TRUE, .N], z[mine %in% TRUE & inkeys %in% TRUE, .N])
R[, np := NP[scen]]
BASE <- R[, c(K6, "est_design", "const", "np"), with = FALSE]
setkeyv(BASE, K6)
ND <- R[, .(nd = sum(est_design), n_const = sum(const), n_rec = .N), by = DPK]
# ---- usability (Section 2.5) and rejections -----------------------------------------------------------------------
# use: fit returned (ok), estimate, SE and p finite, |estimate| < lim, 0 < SE < lim, estimable design, outcome not constant.
# use_conv: the same plus converged == TRUE (an extra condition found in the project code, recorded to measure its effect).
ind <- function(X, lim = 10, wald = TRUE) {
  X <- merge(BASE, X, by = K6, all.x = TRUE)
  if (wald) {
    X[, p_ok := is.finite(estimate) & is.finite(se) & se > 0 & is.finite(estimate / se)]
    X[, use := est_design & !const & ok %in% TRUE & p_ok & abs(estimate) < lim & se < lim & (is.na(p) | is.finite(p))]
    X[, z := estimate / se]
    X[, `:=`(rn = use & 2 * pnorm(-abs(z)) < 0.05, rt = use & 2 * pt(-abs(z), N - np) < 0.05, rt2 = use & 2 * pt(-abs(z), N - 2L) < 0.05)]
  } else {
    X[, use := est_design & !const & ok %in% TRUE & is.finite(p)]
    X[, `:=`(rn = use & p < 0.05, rt = use & p < 0.05, rt2 = use & p < 0.05)]
  }
  X[, use_conv := use & conv %in% TRUE]
  X[, `:=`(rn_conv = use_conv & rn, rt_conv = use_conv & rt)]
  X[, use_nc := est_design & ok %in% TRUE & (if (wald) p_ok & abs(estimate) < lim & se < lim else is.finite(p))]   # without the constant-outcome rule
  X[, rt_nc := use_nc & (if (wald) 2 * pt(-abs(estimate / se), N - np) < 0.05 else p < 0.05)]
  X
}
agg <- function(X, s) X[, .(n = sum(use), kn = sum(rn), kt = sum(rt), kt2 = sum(rt2), n_conv = sum(use_conv), kn_conv = sum(rn_conv), kt_conv = sum(rt_conv),
                            n_nc = sum(use_nc), kt_nc = sum(rt_nc)), by = DPK][, strategy := s]
AGG <- list(); BOOT <- list(); COV <- list()
truth_nat <- fread("tables/V4_S2_truth_glmmnative.csv"); cat("truth file columns:", paste(names(truth_nat), collapse = ","), "\n"); print(head(truth_nat, 3))
add <- function(X, s, boot = FALSE, cov = FALSE, lims = FALSE, wald = TRUE) {
  Y <- ind(X, wald = wald); AGG[[s]] <<- agg(Y, s)
  if (lims) for (L in c(5, 20)) { Z <- ind(X, lim = L, wald = wald); a <- Z[, .(n = sum(use), kt = sum(rt)), by = DPK][, `:=`(strategy = s, lim = L)]; AGG[[paste0(s, "|lim", L)]] <<- a }
  if (boot) BOOT[[s]] <<- Y[, .(dgp, scen, N, ntime, label, sim_id, use, rn, rt)]
  if (cov) COV[[s]] <<- Y[use == TRUE, .(dgp, scen, N, ntime, label, sim_id, estimate, se)]
  invisible(NULL)
}
# ---- study rows (three study mechanisms): mmm (with the check), gee (exch, uncorrected, geepack), gee_md/kc/fg (geesmv), glmer, shrink, glmer_lrt, mmm_nocc
S <- qs_read(file.path(OUT, "x_study.qs2"), nthreads = 2L)
ck("study rows duplicates by key+method: %d", anyDuplicated(S, by = c(K6, "method")))
S[, `:=`(ok = !is.na(term), conv = converged)]
XS <- list()
for (m in c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt", "mmm_nocc")) XS[[paste0("st_", m)]] <- S[method == m, c(K6, "estimate", "se", "p", "ok", "conv"), with = FALSE]
# U0 (stored flag) vs the text definition for the MMM rows
Y <- ind(S[method == "mmm", c(K6, "estimate", "se", "p", "ok", "conv"), with = FALSE]); Y <- merge(Y, S[method == "mmm", c(K6, "U0"), with = FALSE], by = K6, all.x = TRUE)
ck("study mmm: use(text) %d, use(text)&converged %d, stored U0&estimable&!const %d, use(text)&!converged %d",
   Y[, sum(use)], Y[, sum(use_conv)], Y[, sum(U0 %in% TRUE & est_design & !const)], Y[use & !(conv %in% TRUE), .N])
for (m in c("gee", "gee_md", "gee_fg", "glmer", "shrink", "glmer_lrt")) { Y <- ind(S[method == m, c(K6, "estimate", "se", "p", "ok", "conv"), with = FALSE], wald = m != "glmer_lrt")
  ck("study %s: use(text) %d, use&converged %d, use but not converged %d", m, Y[, sum(use)], Y[, sum(use_conv)], Y[use & !(conv %in% TRUE), .N]) }
rm(Y); invisible(gc())
# ---- convergence check reconstructed from its definition (Section 2.4) on the main grid M3-M5 (S6cc), M6-M7 (S6b ns), copula M3-M5 (S7e)
cc_rule <- function(C) { C <- copy(C)
  C[, bnd_mine := A_ok %in% TRUE & (!(A_code %in% c(1L, 2L)) | (is.finite(A_slope_min) & A_slope_min < 1e-3))]
  C[, used_mine := fifelse(bnd_mine & B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se), "B", "A")]   # refit failed (no finite estimate and SE) -> keep A (Section 2.4)
  C[, `:=`(estimate = fifelse(used_mine == "B", B_est, A_est), se = fifelse(used_mine == "B", B_se, A_se), p = fifelse(used_mine == "B", B_p, A_p))]
  C }
C6 <- qs_read(file.path(OUT, "x_s6cc.qs2")); C6 <- C6[, intersect(names(C6), c(K6, "A_ok", "A_has_term", "A_est", "A_se", "A_p", "A_code", "A_slope_min", "bnd", "B_ok", "B_est", "B_se", "B_p", "B_code", "used", "cc_est", "cc_se", "cc_p")), with = FALSE]
if (!"A_has_term" %in% names(C6)) C6[, A_has_term := NA]
C6[is.na(A_has_term), A_has_term := A_ok %in% TRUE & is.finite(A_est)]   # column absent in older S6cc chunks
C7 <- qs_read(file.path(OUT, "x_nscc.qs2")); C8 <- qs_read(file.path(OUT, "x_cop_cc.qs2"))
CCall <- rbind(C6[, src := "main"], C7[, src := "m67"], C8[, src := "copula"], fill = TRUE)
CCall <- cc_rule(CCall)
ck("check: rows %d; stored bnd vs mine disagree %d; stored used vs mine disagree %d; A not ok %d (used B there: %d)", nrow(CCall),
   CCall[bnd != bnd_mine, .N], CCall[used != used_mine, .N], CCall[!(A_ok %in% TRUE), .N], CCall[!(A_ok %in% TRUE) & used == "B", .N])
print(CCall[, .N, by = .(src, A_ok, bnd, B_ok, used)][order(src, -N)])
CCall <- merge(CCall, BASE[, c(K6, "est_design", "const"), with = FALSE], by = K6, all.x = TRUE)
CA <- CCall[, .(n = .N, n_est = sum(est_design %in% TRUE), n_Aok = sum(A_ok %in% TRUE), n_bnd = sum(bnd_mine), n_usedB = sum(used_mine == "B"),
                n_usedB_est = sum(used_mine == "B" & est_design %in% TRUE), n_nlm = sum(A_ok %in% TRUE & !(A_code %in% 1:2)),
                n_slope = sum(A_ok %in% TRUE & A_code %in% 1:2 & is.finite(A_slope_min) & A_slope_min < 1e-3),
                n_nlm_est = sum(est_design %in% TRUE & A_ok %in% TRUE & !(A_code %in% 1:2)), n_bnd_est = sum(est_design %in% TRUE & bnd_mine)), by = c("src", DPK)]
fwrite(CA, file.path(OUT, "cc_agg.csv.gz"))
# stored study mmm (main grid M3-M5 and M6-M7) vs the reconstructed check
z <- merge(S[method == "mmm", c(K6, "estimate", "se"), with = FALSE], CCall[src != "copula", c(K6, "estimate", "se"), with = FALSE], by = K6, suffixes = c("_st", "_cc"))
ck("stored mmm vs reconstructed check (main M3-M5, M6-M7): %d pairs, estimate differs %d, se differs %d", nrow(z),
   z[!(abs(estimate_st - estimate_cc) <= 1e-10 | (is.na(estimate_st) & is.na(estimate_cc))), .N], z[!(abs(se_st - se_cc) <= 1e-8 * pmax(1, abs(se_cc)) | (is.na(se_st) & is.na(se_cc))), .N])
# no-check MMM: fit A on the main grid M3-M5 (S6cc)
XA <- CCall[src == "main", .(dgp, scen, N, ntime, label, sim_id, estimate = A_est, se = A_se, p = A_p, ok = A_ok & A_has_term %in% TRUE, conv = A_code %in% 1:2)]
add(XA, "mmmA_main", cov = FALSE); tm("no-check main")
rm(S, z, XA, C6, C7); invisible(gc())
# ---- copula: study rows + reconstructed check --------------------------------------------------------------------
CS <- qs_read(file.path(OUT, "x_cop_study.qs2")); CS[, `:=`(ok = !is.na(term), conv = converged)]
ck("copula study duplicates: %d", anyDuplicated(CS, by = c(K6, "method")))
for (m in c("gee", "glmer", "shrink", "glmer_lrt")) XS[[paste0("st_", m)]] <- rbind(XS[[paste0("st_", m)]], CS[method == m, c(K6, "estimate", "se", "p", "ok", "conv"), with = FALSE])
Mraw <- CS[method == "mmm", c(K6, "estimate", "se", "p", "ok", "conv"), with = FALSE]
cc3 <- CCall[src == "copula", .(dgp, scen, N, ntime, label, sim_id, e2 = estimate, s2 = se, p2 = p, ok2 = A_ok %in% TRUE & A_has_term %in% TRUE, used_mine)]
Mcc <- merge(Mraw, cc3, by = K6, all = TRUE)
ck("copula mmm: raw rows %d, cc rows %d, raw==A where A used: est differs %d", nrow(Mraw), nrow(cc3), Mcc[used_mine == "A" & ok2 & is.finite(estimate) & abs(estimate - e2) > 1e-10, .N])
Mcc[!is.na(used_mine) & ok2 == TRUE, `:=`(estimate = e2, se = s2, p = p2, ok = TRUE)]
Mcc[!is.na(used_mine) & used_mine == "B", conv := TRUE]
XS[["st_mmm"]] <- rbind(XS[["st_mmm"]], Mcc[, c(K6, "estimate", "se", "p", "ok", "conv"), with = FALSE])
add(Mraw, "mmm_nocheck_cop")
for (s in names(XS)) { m <- sub("^st_", "", s)
  stopifnot(anyDuplicated(XS[[s]], by = K6) == 0L)
  add(XS[[s]], s, boot = m %in% c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt"),
      cov = m %in% c("mmm", "gee_md", "gee_kc", "gee_fg"), lims = m %in% c("mmm", "glmer", "shrink"), wald = m != "glmer_lrt"); tm(s) }
rm(XS); invisible(gc())
rm(CS, Mraw, Mcc, cc3, CCall); invisible(gc()); tm("copula study")
# ---- GEE rows of the S7d and S7e runs (fewclust corrections, CR2) -------------------------------------------------
G <- qs_read(file.path(OUT, "x_gee.qs2"), nthreads = 2L); G[, `:=`(ok = gee_ok, conv = TRUE, p = NA_real_)]
ck("gee rows duplicates: %d", anyDuplicated(G, by = c(K6, "cor")))
for (cn in c("exch", "ind", "ar1")) for (ty in c("unc", "md", "kc", "fg", "cr2")) {
  if (ty == "cr2" && cn == "ar1") next
  X <- G[cor == cn, c(K6, "est", paste0("se_", ty), "ok", "conv", "p"), with = FALSE]; setnames(X, c("est", paste0("se_", ty)), c("estimate", "se"))
  s <- paste0("gee_", cn, "_", ty)
  add(X, s, boot = s %in% c("gee_ar1_fg", "gee_exch_fg", "gee_ind_fg", "gee_exch_kc"), cov = s == "gee_exch_kc", lims = TRUE)
  tm(s)
}
rm(G, X); invisible(gc())
# ---- sign-flip score test ------------------------------------------------------------------------------------------
X <- R[, .(dgp, scen, N, ntime, label, sim_id, estimate = NA_real_, se = NA_real_, p = p_flip, ok = TRUE, conv = TRUE)]
add(X, "flip", boot = TRUE, wald = FALSE); rm(X)
# ---- write ---------------------------------------------------------------------------------------------------------
A <- rbindlist(AGG, fill = TRUE); A <- merge(A, ND, by = DPK); fwrite(A, file.path(OUT, "dp_agg.csv.gz"))
ck("dp_agg rows %d, strategies %d", nrow(A), uniqueN(A$strategy))
qs_save(BOOT, file.path(OUT, "ind_boot.qs2"), nthreads = 2L)
CV <- rbindlist(lapply(names(COV), function(s) COV[[s]][, strategy := s])); qs_save(CV, file.path(OUT, "cov_rows.qs2"), nthreads = 2L)
close(CHK); tm("done")
