
# s10s_rule.R -- S10 "supplemento": independent per-design-point aggregates of every strategy of Table 5 / S26 /
# S28 / S31-S33 / S37 and of the M6-M7 tables (S14, S35), over the 2,376 design points of the four mechanisms.
# Definitions from the rendered text; no project aggregation function is used.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(parallel)}); setDTthreads(3)
RP <- "/home/user/MMM/results_v4"; D <- file.path(RP, "S10/supplemento"); TB <- "/home/user/MMM/tables"
K <- c("dgp","scen","N","ntime","label","sim_id"); CELL <- c("dgp","scen","N","ntime","label"); diag <- list()
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
# ---- universe
R <- rbind(readRDS(file.path(D, "recs_s7d.rds")), readRDS(file.path(D, "recs_s7e.rds")), fill = TRUE)
R[, gen := generated %in% TRUE][, deg := design_degenerate %in% TRUE][, estd := gen & !deg]
R[, const := estd & (sum_y == 0 | sum_y == n_obs)]
diag$universe <- R[, .(dp = uniqueN(paste(dgp, scen, N, ntime)), n = .N, estd = sum(estd), const = sum(const, na.rm = TRUE), notgen_notdeg = sum(!gen & !deg)), by = .(dgp)]
KEYS <- R[estd == TRUE, c(K, "const"), with = FALSE]
NE <- R[, .(n_all = .N, n_est = sum(estd), n_const = sum(const, na.rm = TRUE)), by = CELL]
# ---- M6-M7 MMM fits (A = full latent model, B = refit without time trend) from the S6b chunk files
f <- list.files(file.path(RP, "S6b/ns/rows"), full.names = TRUE)
NSC <- rbindlist(mclapply(f, function(x) readRDS(x)$cc, mc.cores = 3), fill = TRUE); tm("ns cc read")
NS <- readRDS(file.path(RP, "S6b/ns/rows_ns.rds"))
# ---- truths
TT <- fread(file.path(TB, "V4_S4_truth_by_cell.csv"))[, .(dgp, scen, ntime, label, tm = truth_marginal)]
TN <- fread(file.path(TB, "V4_S6b_ns_truth.csv"))[, .(dgp, scen, ntime, label, tm)]
cop_t <- function(scen, label) fifelse(label == "manual", 0, fifelse(scen == "M4", 0.05, 0.2))   # generator (Table S1)
TT <- rbind(TT, TN)
# ---- generic per-strategy aggregation
np_of <- function(scen) fifelse(scen %in% c("M2","M5"), 5L, 4L)
agg <- function(W, strat, marg = TRUE) {     # W: K + estimate, se, ok (fit returned), [p for lrt/flip]
  W <- merge(W, KEYS, by = K)               # estimable designs only
  W[, np := np_of(scen)]
  W[, fin10 := ok %in% TRUE & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10]
  W[, fin5 := fin10 & abs(estimate) < 5 & se < 5]
  W[, fin20 := ok %in% TRUE & is.finite(estimate) & is.finite(se) & abs(estimate) < 20 & se > 0 & se < 20]
  W[, z := estimate / se]
  W[, `:=`(rt = abs(z) > qt(0.975, N - np), rt2 = abs(z) > qt(0.975, N - 2), rn = abs(z) > qnorm(0.975))]
  if (marg) {
    W <- merge(W, TT, by = c("dgp","scen","ntime","label"), all.x = TRUE)
    W[dgp == "Copula-AR1", tm := cop_t(scen, label)]
    W[, `:=`(cvn = abs(estimate - tm) <= qnorm(0.975) * se, cvt = abs(estimate - tm) <= qt(0.975, N - np) * se)]
  } else W[, `:=`(tm = NA_real_, cvn = NA, cvt = NA)]
  W[, use := fin10 & !const]
  A <- W[, .(n_row = .N, n_use = sum(use), k_t = sum(use & rt), k_t2 = sum(use & rt2), k_n = sum(use & rn),
             n_use5 = sum(fin5 & !const), k_t5 = sum(fin5 & !const & rt), n_use20 = sum(fin20 & !const), k_t20 = sum(fin20 & !const & rt),
             n_inc = sum(fin10), k_inc = sum(fin10 & rt),
             s_est = sum(estimate[use]), s_est2 = sum(estimate[use]^2), s_se = sum(se[use]), s_err = sum((estimate - tm)[use]),
             cov_n = sum(use & cvn), cov_t = sum(use & cvt), tm = tm[1]), by = CELL]
  A[, strategy := strat]; A
}
agg_p <- function(W, strat) {               # p-value tests (likelihood ratio, sign-flip)
  W <- merge(W, KEYS, by = K); W[, use := ok %in% TRUE & is.finite(p) & !const]; W[, inc := ok %in% TRUE & is.finite(p)]
  A <- W[, .(n_row = .N, n_use = sum(use), k_t = sum(use & p < 0.05), k_t_le = sum(use & p <= 0.05), n_inc = sum(inc), k_inc = sum(inc & p < 0.05)), by = CELL]
  A[, strategy := strat]; A
}
OUT <- list()
# ---- GEE strategies (fewclust SEs, three working correlations, both GEE row sources)
G <- rbind(readRDS(file.path(D, "gee_s7d.rds")), readRDS(file.path(D, "gee_s7e.rds")), fill = TRUE); tm("gee loaded")
diag$smv_vs_fc <- G[cor == "exch", .(n_md = sum(is.finite(smv_md) & is.finite(se_md)), md_1pct = sum(abs(smv_md / se_md - 1) < 0.01, na.rm = TRUE),
                                    n_fg = sum(is.finite(smv_fg) & is.finite(se_fg)), fg_1pct = sum(abs(smv_fg / se_fg - 1) < 0.01, na.rm = TRUE),
                                    n_kc = sum(is.finite(smv_kc) & is.finite(se_kc)), kc_1pct = sum(abs(smv_kc / se_kc - 1) < 0.01, na.rm = TRUE)), by = .(src = dgp == "Copula-AR1")]
G[, se_avg := sqrt((se_md^2 + se_kc^2) / 2)]
for (cr in c("exch","ind","ar1")) for (co in c("unc","md","kc","fg","cr2","avg")) {
  if (co == "cr2" && cr == "ar1") next
  W <- G[cor == cr, c(K, "est", paste0("se_", co), "gee_ok"), with = FALSE]; setnames(W, c(K, "estimate", "se", "ok"))
  OUT[[paste0("gee_", cr, "_", co)]] <- agg(W, paste0("gee_", cr, "_", co))
}
# geesmv versions (exchangeable) for Table S14 rows and the Section S4 comparison
for (co in c("md","fg","kc")) { W <- G[cor == "exch", c(K, "est", paste0("smv_", co), "gee_ok"), with = FALSE]; setnames(W, c(K, "estimate", "se", "ok"))
  OUT[[paste0("gee_exch_smv_", co)]] <- agg(W, paste0("gee_exch_smv_", co)) }
rm(G, W); gc(); tm("gee strategies done")
# ---- MMM with the check, rebuilt from the text definition, and without the check
cc_rule <- function(C) { C[, bnd_text := !(A_code %in% c(1L, 2L)) | (is.finite(A_slope_min) & A_slope_min < 0.001)]
  C[, useB := bnd_text & (B_ok %in% TRUE)]
  C[, .(dgp, scen, N, ntime, label, sim_id, estimate = fifelse(useB, B_est, A_est), se = fifelse(useB, B_se, A_se),
        ok = fifelse(useB, TRUE, A_ok %in% TRUE), useB, bnd_text, used)] }
S <- readRDS(file.path(D, "study_cc.rds"))
CCm <- cc_rule(readRDS(file.path(RP, "S6cc/cc_all.rds")))
E <- readRDS(file.path(D, "study_s7e.rds")); CE <- readRDS(file.path(D, "cc_s7e.rds"))
diag$cop_cc_cols <- names(CE)
CCe <- cc_rule(CE); CCn <- cc_rule(NSC)
diag$cc_used_agree <- rbind(CCm[, .(src = "S6cc", n = .N, disagree = sum((useB) != (used == "B"), na.rm = TRUE))],
                            CCe[, .(src = "S7e", n = .N, disagree = sum((useB) != (used == "B"), na.rm = TRUE))],
                            CCn[, .(src = "S6b", n = .N, disagree = sum((useB) != (used == "B"), na.rm = TRUE))])
diag$ns_intervened <- CCn[, .(n = .N, bnd = sum(bnd_text)), by = .(dgp, scen, label)]
mm_simple <- function(X) X[, .(dgp, scen, N, ntime, label, sim_id, estimate, se, ok = converged %in% TRUE & is.finite(p))]
MMc <- rbind(mm_simple(S[method == "mmm" & scen %in% c("M1","M2","M1c")]), CCm[, c(K, "estimate","se","ok"), with = FALSE],
             mm_simple(E[method == "mmm" & scen %in% c("M1","M2")]), CCe[, c(K, "estimate","se","ok"), with = FALSE],
             CCn[, c(K, "estimate","se","ok"), with = FALSE])
OUT$mmm_cc <- agg(MMc, "mmm_cc")
MMn <- rbind(mm_simple(readRDS(file.path(D, "study_mmm_nocc.rds"))), mm_simple(E[method == "mmm"]), mm_simple(NS[method == "mmm_nocc"]))
OUT$mmm_nocc <- agg(MMn, "mmm_nocc"); rm(MMc, MMn); gc()
# stored check (rows_ns mmm) for M6-M7 as a cross-check of the rebuilt check
OUT$mmm_cc_storedns <- agg(mm_simple(NS[method == "mmm"]), "mmm_cc_storedns")
# ---- mixed model, shrinkage, likelihood ratio, geesmv corrections from the study rows
for (m in c("glmer","shrink","gee_md","gee_fg","gee")) {
  W <- rbind(mm_simple(S[method == m]), mm_simple(E[method == m]), mm_simple(NS[method == m]))
  OUT[[paste0("study_", m)]] <- agg(W, paste0("study_", m), marg = !(m %in% c("glmer"))) }
W <- rbind(S[method == "glmer_lrt", .(dgp, scen, N, ntime, label, sim_id, p, ok = converged)], E[method == "glmer_lrt", .(dgp, scen, N, ntime, label, sim_id, p, ok = converged)],
           NS[method == "glmer_lrt", .(dgp, scen, N, ntime, label, sim_id, p, ok = converged)])
OUT$glmer_lrt <- agg_p(W, "glmer_lrt")
# GLMM SE ratio on Copula (text S8): per design point model SE / SD of estimates, all usable
W <- E[method == "glmer"]; W <- merge(W, KEYS, by = K)
W[, use := converged %in% TRUE & is.finite(estimate) & is.finite(se) & is.finite(p) & abs(estimate) < 10 & se > 0 & se < 10 & !const]
CG <- W[use == TRUE, .(n = .N, mse = mean(se), sde = sd(estimate), mest = mean(estimate)), by = CELL]
fwrite(CG, file.path(D, "s10s_cop_glmm_ser.csv"))
rm(S, E, NS, W); gc()
# ---- sign-flip test
W <- R[, .(dgp, scen, N, ntime, label, sim_id, p = p_flip, ok = TRUE)]; OUT$flip <- agg_p(W, "flip")
W <- R[, .(dgp, scen, N, ntime, label, sim_id, p = p_flip_strict, ok = TRUE)]; OUT$flip_strict <- agg_p(W, "flip_strict")
A <- rbindlist(OUT, fill = TRUE); A <- merge(A, NE, by = CELL, all.x = TRUE)
fwrite(A, file.path(D, "s10s_rule_cells.csv"))
# ---- M6-M7 check intervention per cell, with and without (Table S35)
X <- merge(CCn, KEYS, by = K); X <- merge(X, TN, by = c("dgp","scen","ntime","label"), all.x = TRUE)
X[, use := ok & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & !const]
X[, cvt := abs(estimate - tm) <= qt(0.975, N - 4L) * se]
B35 <- X[, .(n_use = sum(use), n_bnd = sum(use & useB), n_bnd_any = sum(use & bnd_text), cov_with = sum(use & useB & cvt), cov_wo = sum(use & !useB & cvt),
             s_se_with = sum(se[use & useB]), s_se_wo = sum(se[use & !useB]), sd_all = sd(estimate[use]), sd_with = sd(estimate[use & useB]), sd_wo = sd(estimate[use & !useB]),
             s_err = sum((estimate - tm)[use]), cov_t = sum(use & cvt), s_se = sum(se[use]), tm = tm[1]), by = CELL]
fwrite(B35, file.path(D, "s10s_ns_cc_cells.csv"))
sink(file.path(D, "s10s_rule_diag.txt")); print(diag); sink()
tm("RULE OK")
