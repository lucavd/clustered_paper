# s10t_ccvar.R -- S10 "Testo" follow-up: three variants of the convergence check of the marginalized model.
# (i) as today: at the boundary, refit B used when it has finite estimate and SE (any nlm code), else fit A;
# (ii) B used only if nlm reports convergence for B (code 1 or 2) and estimate/SE finite, else fit A;
# (iii) B used only if it converges (code 1 or 2, finite); otherwise no usable result.
# (ii3), (iii3): as (ii), (iii) with nlm code 1, 2 or 3 accepted (3 = approximate minimum).
# Boundary rule (Section 2.4) recomputed from A_code and A_slope_min. Scenarios without the check (M1, M2, M1c) keep the
# stored study rows. Usability and rejections as in s10t_dp.R (Section 2.5 + constant-outcome rule). Read-only inputs:
# results_v4/S10/testo/x_*.qs2 (copies of S6cc, S6b-ns cc, S7e cc/study, study rows, replicate records).
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L); setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); DPK <- c("dgp", "scen", "N", "ntime", "label")
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
R <- qs_read(file.path(OUT, "x_recs.qs2")); R[, est_design := generated & !design_degenerate]; R[, const := est_design & (sum_y == 0 | sum_y == n_obs)]
BASE <- R[, c(K6, "est_design", "const"), with = FALSE]
# stored rows for the scenarios without the check
S <- qs_read(file.path(OUT, "x_study.qs2"))[method == "mmm" & scen %in% c("M1", "M2", "M1c"), c(K6, "term", "estimate", "se", "p"), with = FALSE]
CS <- qs_read(file.path(OUT, "x_cop_study.qs2"))[method == "mmm" & scen %in% c("M1", "M2"), c(K6, "term", "estimate", "se", "p"), with = FALSE]
NOCK <- rbind(S, CS)[, ok := !is.na(term)][, term := NULL]
cols <- c(K6, "A_ok", "A_est", "A_se", "A_p", "A_code", "A_slope_min", "B_ok", "B_est", "B_se", "B_p", "B_code", "used")
C6 <- qs_read(file.path(OUT, "x_s6cc.qs2")); C7 <- qs_read(file.path(OUT, "x_nscc.qs2")); C8 <- qs_read(file.path(OUT, "x_cop_cc.qs2"))
CC <- rbind(C6[, ..cols][, src := "main"], C7[, ..cols][, src := "m67"], C8[, ..cols][, src := "copula"])
CC[, A_fin := A_ok %in% TRUE & is.finite(A_est)]
CC[, bnd := A_ok %in% TRUE & (!(A_code %in% 1:2) | (is.finite(A_slope_min) & A_slope_min < 1e-3))]
CC[, B_fin := B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se)]
CC[, B_conv := B_fin & B_code %in% 1:2]; CC[, B_conv3 := B_fin & B_code %in% 1:3]
cat("rows", nrow(CC), "; at boundary", CC[bnd == TRUE, .N], "; B finite", CC[bnd & B_fin, .N], "; B finite and converged", CC[bnd & B_conv, .N], "\n")
print(CC[bnd == TRUE, .(n = .N, B_fin = sum(B_fin), B_conv = sum(B_conv), B_fin_notconv = sum(B_fin & !B_conv)), by = src])
pick <- function(v) {
  X <- copy(CC)
  if (v == "i") X[, choice := fifelse(!bnd, "A", fifelse(B_fin, "B", "A"))]
  if (v == "ii") X[, choice := fifelse(!bnd, "A", fifelse(B_conv, "B", "A"))]
  if (v == "iii") X[, choice := fifelse(!bnd, "A", fifelse(B_conv, "B", "none"))]
  if (v == "ii3") X[, choice := fifelse(!bnd, "A", fifelse(B_conv3, "B", "A"))]
  if (v == "iii3") X[, choice := fifelse(!bnd, "A", fifelse(B_conv3, "B", "none"))]
  X[, `:=`(estimate = fifelse(choice == "B", B_est, fifelse(choice == "A", A_est, NA_real_)),
           se = fifelse(choice == "B", B_se, fifelse(choice == "A", A_se, NA_real_)),
           p = fifelse(choice == "B", B_p, fifelse(choice == "A", A_p, NA_real_)),
           ok = fifelse(choice == "B", TRUE, fifelse(choice == "A", A_fin, FALSE)))]
  if (v == "i") cat("variant (i) choice vs stored 'used' disagree:", X[!is.na(used) & choice != used, .N], "of", X[!is.na(used), .N], "\n")
  rbind(NOCK, X[, c(K6, "estimate", "se", "p", "ok"), with = FALSE])
}
OUTL <- list()
for (v in c("i", "ii", "iii", "ii3", "iii3")) {
  X <- merge(BASE, pick(v), by = K6, all.x = TRUE); stopifnot(nrow(X) == nrow(BASE))
  X[, np := NP[scen]]
  X[, use := est_design & !const & ok %in% TRUE & is.finite(estimate) & is.finite(se) & se > 0 & abs(estimate) < 10 & se < 10 & (is.na(p) | is.finite(p))]
  X[, z := estimate / se]
  X[, `:=`(rn = use & 2 * pnorm(-abs(z)) < 0.05, rt = use & 2 * pt(-abs(z), N - np) < 0.05)]
  OUTL[[v]] <- X[, .(n = sum(use), kn = sum(rn), kt = sum(rt), nd = sum(est_design)), by = DPK][, variant := v]
  cat("variant", v, ": usable", X[, sum(use)], "\n")
}
fwrite(rbindlist(OUTL), file.path(OUT, "ccvar_dp.csv"))
