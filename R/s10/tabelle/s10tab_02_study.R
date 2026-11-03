
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R"); setDTthreads(2L)
# Script 2: study rows of the main grid and M1c (variant_cc row set, before the D-82 and kcfc edits), convergence
# check rebuilt from results_v4/S6cc/cc_all.rds, KC (fewclust) from script 1. Universe = S7d recs (script 1).
t0 <- Sys.time()
R <- readRDS(file.path(OUT, "s1_recs.rds"))[!(scen %in% c("M6", "M7"))]
cat("universe", nrow(R), "estimable", sum(R$estim), "const", sum(R$const), "\n")
M <- c("mmm", "gee", "gee_md", "gee_fg", "glmer", "shrink", "glmer_lrt")
V <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))
V <- V[method %in% M & term %in% "time:trt", .(dgp, scen, N, ntime, label, sim_id, method, estimate, se, p, stat, converged)]
invisible(gc()); V[, scen := norm_scen(scen)]
cat(sprintf("rows kept %d (%.1f min)\n", nrow(V), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
dup <- V[, .N, by = c(K6, "method")][N > 1L]; cat("duplicated (key, method):", nrow(dup), "\n")
# consistency of stored p with the stored statistic (information only)
chk <- V[, .(n = .N, wald_p_mismatch = sum(method != "glmer_lrt" & is.finite(p) & is.finite(estimate / se) & abs(p - 2 * pnorm(-abs(estimate / se))) > 1e-6),
             lrt_p_mismatch = sum(method == "glmer_lrt" & is.finite(p) & is.finite(stat) & abs(p - pchisq(stat, 1, lower.tail = FALSE)) > 1e-6)), by = method]
print(chk)
CC <- as.data.table(readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))); CC[, scen := norm_scen(scen)]
O <- cbind(CC[, .(dgp, scen, N, ntime, label)], own_cc(CC)); rm(CC)
print(O[, .(n = .N, bnd_disagree = sum(xor(bnd_own %in% TRUE, bnd_proj %in% TRUE)), used_disagree = sum(used_own != used_proj, na.rm = TRUE), used_proj_na = sum(is.na(used_proj))), by = scen])
# project's mmm row on M3-M5 versus the rebuilt check (information: which reading the stored row follows)
cmpm <- merge(V[method == "mmm" & scen %in% c("M3", "M4", "M5"), c(K6, "estimate", "se"), with = FALSE][, .SD, .SDcols = c(K6, "estimate", "se")][, se_row := se][, se := NULL], O, by = K6, all = TRUE)
print(cmpm[, .(n = .N, eq_literal = sum(abs(estimate - est) < 1e-10 & abs(se_row - se) < 1e-10, na.rm = TRUE),
               eq_finite = sum(abs(estimate - est_f) < 1e-10 & abs(se_row - se_f) < 1e-10, na.rm = TRUE), row_missing = sum(is.na(estimate)),
               differ_where_used_disagree = sum(used_own != used_proj & !(abs(estimate - est) < 1e-10), na.rm = TRUE))])
mk <- function(X, s, lrt = FALSE) {
  X <- merge(R[, c(K6, "estim", "const"), with = FALSE], X, by = K6, all.x = TRUE); base <- X$estim & !X$const
  if (lrt) { use <- base & (X$ok %in% TRUE) & is.finite(X$p); rn <- rt <- use & X$p < 0.05 } else {
    use <- base & wald_ok(X$ok, X$est, X$se); rr <- zrej(use, X$est, X$se, X$N - P_OF[X$scen]); rn <- rr$rn %in% TRUE; rt <- rr$rt %in% TRUE }
  data.table(X[, K6, with = FALSE], strategy = s, use = use, rn = rn, rt = rt) }
Vs <- function(m) V[method == m, c(K6, "estimate", "se", "p", "converged"), with = FALSE][, .SD, .SDcols = c(K6, "estimate", "se", "p", "converged")]
ren <- function(X) setnames(X, c("estimate", "converged"), c("est", "ok"))
IND <- list()
for (m in setdiff(M, "mmm")) IND[[m]] <- mk(ren(Vs(m)), m, lrt = (m == "glmer_lrt"))
IND[["mmm_projrow"]] <- mk(ren(Vs("mmm")), "mmm_projrow")
mm12 <- ren(Vs("mmm"))[!(scen %in% c("M3", "M4", "M5"))]
IND[["mmm_cc"]] <- mk(rbind(mm12, O[, .(dgp, scen, N, ntime, label, sim_id, est, se, p = NA_real_, ok)]), "mmm_cc")
IND[["mmm_cc_f"]] <- mk(rbind(mm12, O[, .(dgp, scen, N, ntime, label, sim_id, est = est_f, se = se_f, p = NA_real_, ok = ok_f)]), "mmm_cc_f")
KC <- readRDS(file.path(OUT, "s1_kc_rep.rds")); KC[, rt := FALSE]
IND[["gee_exch_kc_n"]] <- merge(R[scen != "M1c", K6, with = FALSE], KC, by = K6, all.x = TRUE)[, strategy := "gee_exch_kc"][is.na(use), `:=`(use = FALSE, rn = FALSE)]
rm(V); invisible(gc())
I <- rbindlist(IND, use.names = TRUE); rm(IND); invisible(gc())
ND <- R[, .(nd = sum(estim), nconst = sum(const), nrecs = .N), by = KD]
DP <- merge(I[, .(n = sum(use), kn = sum(rn), kt = sum(rt)), by = c(KD, "strategy")], ND, by = KD)
saveRDS(DP, file.path(OUT, "s2_study_dp.rds"))
B <- I[label == "manual" & scen %in% c("M1", "M2", "M3", "M4", "M5"), .(u = sum(use), r = sum(rn)), by = .(dgp, scen, N, sim_id, strategy)]
saveRDS(B, file.path(OUT, "s2_boot_rep.rds"))
cat(sprintf("done %.1f min; dp %d rows; boot %d rows\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(DP), nrow(B)))
