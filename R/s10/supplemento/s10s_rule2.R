
# s10s_rule2.R -- S10 "supplemento", follow-up: (1) MMM with the check as stored by the project (to measure the effect
# of the divergence between the stored check and the text definition), (2) sign-flip test with the rules of main 2.4
# (p0_guard), (3) width ratio of the inverted sign-flip interval, (4) AGQ row counts, (5) M6-M7 per-point measures
# for Table S35 under both versions of the check.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(parallel)}); setDTthreads(2)
RP <- "/home/user/MMM/results_v4"; D <- file.path(RP, "S10/supplemento"); TB <- "/home/user/MMM/tables"
K <- c("dgp","scen","N","ntime","label","sim_id"); CELL <- c("dgp","scen","N","ntime","label")
R <- rbind(readRDS(file.path(D, "recs_s7d.rds")), readRDS(file.path(D, "recs_s7e.rds")), fill = TRUE)
R[, estd := generated %in% TRUE & !(design_degenerate %in% TRUE)][, const := estd & (sum_y == 0 | sum_y == n_obs)]
KEYS <- R[estd == TRUE, c(K, "const"), with = FALSE]
np_of <- function(scen) fifelse(scen %in% c("M2","M5"), 5L, 4L)
cop_t <- function(scen, label) fifelse(label == "manual", 0, fifelse(scen == "M4", 0.05, 0.2))
TT <- rbind(fread(file.path(TB, "V4_S4_truth_by_cell.csv"))[, .(dgp, scen, ntime, label, tm = truth_marginal)],
            fread(file.path(TB, "V4_S6b_ns_truth.csv"))[, .(dgp, scen, ntime, label, tm)])
agg <- function(W) { W <- merge(W, KEYS, by = K); W[, np := np_of(scen)]
  W <- merge(W, TT, by = c("dgp","scen","ntime","label"), all.x = TRUE); W[dgp == "Copula-AR1", tm := cop_t(scen, label)]
  W[, fin := ok %in% TRUE & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10]
  W[, fin5 := fin & abs(estimate) < 5 & se < 5]; W[, fin20 := ok %in% TRUE & is.finite(estimate) & is.finite(se) & abs(estimate) < 20 & se > 0 & se < 20]
  W[, use := fin & !const]; W[, z := estimate / se]
  W[, `:=`(rt = abs(z) > qt(0.975, N - np), rt2 = abs(z) > qt(0.975, N - 2), rn = abs(z) > qnorm(0.975))]
  W[, `:=`(cvn = abs(estimate - tm) <= qnorm(0.975) * se, cvt = abs(estimate - tm) <= qt(0.975, N - np) * se)]
  W[, .(n_use = sum(use), k_t = sum(use & rt), k_t2 = sum(use & rt2), k_n = sum(use & rn),
        n_use5 = sum(fin5 & !const), k_t5 = sum(fin5 & !const & rt), n_use20 = sum(fin20 & !const), k_t20 = sum(fin20 & !const & rt),
        s_est = sum(estimate[use]), s_est2 = sum(estimate[use]^2), s_se = sum(se[use]), s_se2 = sum(se[use]^2), s_err = sum((estimate - tm)[use]),
        cov_n = sum(use & cvn), cov_t = sum(use & cvt)), by = CELL] }
# (1) stored check
S <- readRDS(file.path(D, "study_cc.rds"))[method == "mmm", .(dgp, scen, N, ntime, label, sim_id, estimate, se, ok = converged %in% TRUE & is.finite(p))]
NS <- readRDS(file.path(RP, "S6b/ns/rows_ns.rds"))[method == "mmm", .(dgp, scen, N, ntime, label, sim_id, estimate, se, ok = converged %in% TRUE & is.finite(p))]
E <- readRDS(file.path(D, "study_s7e.rds")); CE <- readRDS(file.path(D, "cc_s7e.rds"))
E1 <- E[method == "mmm" & scen %in% c("M1","M2"), .(dgp, scen, N, ntime, label, sim_id, estimate, se, ok = converged %in% TRUE & is.finite(p))]
E2 <- CE[, .(dgp, scen, N, ntime, label, sim_id, estimate = cc_est, se = cc_se, ok = is.finite(cc_p))]
A1 <- agg(rbind(S, NS, E1, E2)); A1[, strategy := "mmm_cc_stored"]
# Copula MMM-with-check per the text, from the S7e A/B fits, also through this function (for S28 RMS SE)
CE[, bnd_text := !(A_code %in% c(1L, 2L)) | (is.finite(A_slope_min) & A_slope_min < 0.001)][, useB := bnd_text & (B_ok %in% TRUE)]
E3 <- CE[, .(dgp, scen, N, ntime, label, sim_id, estimate = fifelse(useB, B_est, A_est), se = fifelse(useB, B_se, A_se), ok = fifelse(useB, TRUE, A_ok %in% TRUE))]
A1b <- agg(rbind(E1, E3)); A1b[, strategy := "mmm_cc_text_cop"]
rm(S, E, CE); gc()
# (2) sign-flip with the main 2.4 rules (p0_guard), all 2,376 design points
C <- readRDS(file.path(D, "ci_s7d2.rds")); C <- merge(C[, c(K, "p0_guard", "lo", "hi", "ar1_gee_ok", "ar1_est", "ar1_se_fg"), with = FALSE], KEYS, by = K)
A2 <- C[, .(n_use = sum(is.finite(p0_guard) & !const), k_t = sum(is.finite(p0_guard) & !const & p0_guard < 0.05),
            k_le = sum(is.finite(p0_guard) & !const & p0_guard <= 0.05), n_inc = sum(is.finite(p0_guard))), by = CELL]; A2[, strategy := "flip_guard"]
NE <- R[, .(n_all = .N, n_est = sum(estd)), by = CELL]
fwrite(merge(rbind(A1, A1b, A2, fill = TRUE), NE, by = CELL, all.x = TRUE), file.path(D, "s10s_rule2_cells.csv"))
# (3) width ratio of the inverted interval
C[, np := np_of(scen)]
C[, gu := ar1_gee_ok %in% TRUE & is.finite(ar1_est) & is.finite(ar1_se_fg) & abs(ar1_est) < 10 & ar1_se_fg > 0 & ar1_se_fg < 10 & !const]
C[, `:=`(fw = hi - lo, gw = 2 * qt(0.975, N - np) * ar1_se_fg, gwn = 2 * qnorm(0.975) * ar1_se_fg)]
B <- C[is.finite(lo) & is.finite(hi) & gu == TRUE]
gm <- function(r) 100 * (exp(mean(log(r))) - 1)
W3 <- rbindlist(list(
  B[, .(set = "all", n = .N, n_fw0 = sum(fw <= 0), gm_t = gm(fw[fw > 0] / gw[fw > 0]), gm_n = gm(fw[fw > 0] / gwn[fw > 0])), by = label],
  B[N >= 75, .(set = "N>=75", n = .N, n_fw0 = sum(fw <= 0), gm_t = gm(fw[fw > 0] / gw[fw > 0]), gm_n = gm(fw[fw > 0] / gwn[fw > 0])), by = label]))
W3dp <- B[fw > 0, .(n = .N, lr = mean(log(fw / gw))), by = CELL]
fwrite(W3, file.path(D, "s10s_ciwidth.csv")); fwrite(W3dp, file.path(D, "s10s_ciwidth_dp.csv"))
# finite share in small designs, pooled vs mean over design points (Section S7, 14%)
Cs <- C[label == "mmlb" & ntime <= 3 & N <= 25]
fw2 <- data.table(item = c("pooled", "dpmean", "pooled_nc"), value = c(100 * Cs[, mean(is.finite(lo) & is.finite(hi))],
        100 * Cs[, .(f = mean(is.finite(lo) & is.finite(hi))), by = CELL][, mean(f)], 100 * Cs[const == FALSE, mean(is.finite(lo) & is.finite(hi))]))
fwrite(fw2, file.path(D, "s10s_ci_small.csv")); rm(C, B); gc()
# (4) AGQ counts
Wg <- readRDS(file.path(D, "rows_s7g.rds"))
fwrite(Wg[, .(n = .N, conv = sum(converged %in% TRUE), fin = sum(is.finite(estimate))), by = nAGQ], file.path(D, "s10s_agq_counts.csv"))
# (5) M6-M7 per point (Table S35): both versions of the check
f <- list.files(file.path(RP, "S6b/ns/rows"), full.names = TRUE)
NSC <- rbindlist(mclapply(f, function(x) readRDS(x)$cc, mc.cores = 2), fill = TRUE)
NSC[, bnd_text := !(A_code %in% c(1L, 2L)) | (is.finite(A_slope_min) & A_slope_min < 0.001)][, useB := bnd_text & (B_ok %in% TRUE)]
TN <- fread(file.path(TB, "V4_S6b_ns_truth.csv"))[, .(dgp, scen, ntime, label, tm)]
out <- list()
for (v in c("text", "stored")) {
  X <- copy(NSC)
  if (v == "text") X[, `:=`(e = fifelse(useB, B_est, A_est), s = fifelse(useB, B_se, A_se), ok = fifelse(useB, TRUE, A_ok %in% TRUE), w = useB)]
  else X[, `:=`(e = cc_est, s = cc_se, ok = is.finite(cc_p), w = used == "B")]
  X <- merge(X, KEYS, by = K); X <- merge(X, TN, by = c("dgp","scen","ntime","label"), all.x = TRUE)
  X[, use := ok %in% TRUE & is.finite(e) & is.finite(s) & abs(e) < 10 & s > 0 & s < 10 & !const]
  X[, cvt := abs(e - tm) <= qt(0.975, N - 4L) * s]
  P <- X[use == TRUE, .(n = .N, bias = mean(e - tm), sde = sd(e), rms = sqrt(mean(s^2)), mse = mean(s), cov = 100 * mean(cvt),
                        n_w = sum(w %in% TRUE), cov_w = 100 * mean(cvt[w %in% TRUE]), cov_wo = 100 * mean(cvt[!(w %in% TRUE)]),
                        rms_w = sqrt(mean(s[w %in% TRUE]^2)), rms_wo = sqrt(mean(s[!(w %in% TRUE)]^2)),
                        sde_w = sd(e[w %in% TRUE]), sde_wo = sd(e[!(w %in% TRUE)]), k_nc = sum(!cvt)), by = CELL]
  P[, version := v]; out[[v]] <- P }
fwrite(rbindlist(out), file.path(D, "s10s_ns_s35_dp.csv"))
cat("RULE2 OK\n")
