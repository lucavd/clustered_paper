# Provenance (2026-11-02): read-only check run from a scratch directory (writes to ./out); see docs/steps/S9.md and docs/decisions/D-82.md.

# const_check2.R -- read-only: constant-outcome datasets among the usable results of the marginalized model with the
# check (mmm_cc) and of the sign-flip test with the D-57 rules; Table 5 quantities with and without them; gate vs Table 5.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L)
REPO <- "/home/user/MMM"; dir.create("out", showWarnings = FALSE)
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
RC0 <- c("sim_id", "generated", "design_degenerate", "n_obs", "sum_y", "dgp", "scen", "N", "ntime", "label")
fd <- list.files(file.path(REPO, "results_v4/S7d/rows"), full.names = TRUE, pattern = "\\.rds$")
fe <- list.files(file.path(REPO, "results_v4/S7e/rows"), full.names = TRUE, pattern = "\\.rds$")
RECS <- rbindlist(parallel::mclapply(fd, function(f) as.data.table(readRDS(f)$recs)[, ..RC0], mc.cores = 16L))
LE <- parallel::mclapply(fe, function(f) { x <- readRDS(f)
  st <- as.data.table(x$study)[method == "mmm" & term %in% "time:trt", .(dgp, scen, N, ntime, label, sim_id, estimate, se, p, converged)]
  cc <- as.data.table(x$cc); if (nrow(cc)) cc <- cc[A_ok == TRUE & A_has_term == TRUE, .(dgp, scen, N, ntime, label, sim_id, cc_est, cc_se, cc_p)]
  list(recs = as.data.table(x$recs)[, ..RC0], st = st, cc = cc) }, mc.cores = 16L)
RECE <- rbindlist(lapply(LE, `[[`, "recs")); STE <- rbindlist(lapply(LE, `[[`, "st")); CCE <- rbindlist(lapply(LE, `[[`, "cc"), fill = TRUE); rm(LE)
cat("copula dgp labels:", unique(RECE$dgp), "| S7d dgp:", unique(RECS$dgp), "\n")
RECA <- rbind(RECS, RECE); RECA[, const := generated & !design_degenerate & (sum_y == 0 | sum_y == n_obs)]
cat(sprintf("recs %d (unique %d), constant %d\n", nrow(RECA), uniqueN(RECA[, ..K6]), RECA[, sum(const)]))
# ---- marginalized model with the check
u0 <- function(converged, estimate, se, p) { w <- converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p); w[is.na(w)] <- FALSE; w }
KEEP <- c(K6, "method", "term", "estimate", "se", "U0")
SA <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))[method == "mmm" & term %in% "time:trt", ..KEEP]
SB <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds"))[method == "mmm" & term %in% "time:trt", ..KEEP]
M1 <- rbind(SA, SB)[, .(dgp, scen, N, ntime, label, sim_id, estimate, se, U0)]; rm(SA, SB)
Sm <- merge(STE, CCE, by = K6, all.x = TRUE)
Sm[!is.na(cc_est), `:=`(estimate = cc_est, se = cc_se, p = cc_p)]
Sm[, U0 := u0(converged, estimate, se, p)]
M <- rbind(M1, Sm[, .(dgp, scen, N, ntime, label, sim_id, estimate, se, U0)])
M <- merge(RECA[, c(K6, "design_degenerate", "const"), with = FALSE], M, by = K6, all.x = TRUE)
M[, np := NP[scen]]; M[, use := !design_degenerate & U0 %in% TRUE]; M[, rej := use & 2 * pt(-abs(estimate / se), N - np) < 0.05]
# ---- sign-flip test with the D-57 rules
source(path.expand("~/MMM/R/s7d2/s7d2_common.R")); TA <- s7d2_tasks()
FX <- rbindlist(parallel::mclapply(TA$file, function(f) readRDS(f)[, .(dgp, scen, N, ntime, label, sim_id, p0_guard)], mc.cores = 16L))
cat("flip dgp labels:", unique(FX$dgp), "\n")
FX <- merge(RECA[, c(K6, "design_degenerate", "const"), with = FALSE], FX, by = K6, all.x = TRUE)
FX[, use := !design_degenerate & is.finite(p0_guard)]; FX[, rej := use & p0_guard < 0.05]
IND <- list(mmm_cc = M[, c(K6, "use", "rej", "const"), with = FALSE], flip = FX[, c(K6, "use", "rej", "const"), with = FALSE])
ND <- RECA[, .(nd = sum(!design_degenerate)), by = .(dgp, scen, N, ntime, label)]
cp_lo <- function(k, n) fifelse(k > 0, qbeta(0.025, k, n - k + 1), 0); cp_hi <- function(k, n) fifelse(k < n, qbeta(0.975, k + 1, n - k), 1)
DP <- rbindlist(lapply(names(IND), function(s) IND[[s]][, .(n = sum(use), k = sum(rej), n_c = sum(use & const), k_c = sum(rej & const)), by = .(dgp, scen, N, ntime, label)][, strategy := s]))
DP <- merge(DP, ND, by = c("dgp", "scen", "N", "ntime", "label"))
H0 <- DP[label == "manual"]; H1 <- DP[label == "mmlb"]
for (v in c("orig", "corr")) { nn <- if (v == "orig") H0$n else H0$n - H0$n_c; kk <- if (v == "orig") H0$k else H0$k - H0$k_c
  H0[, (paste0("pu_", v)) := 100 * nn / nd]; H0[, (paste0("hi_", v)) := nn >= 20 & cp_lo(kk, nn) > 0.075]
  H0[, (paste0("lo_", v)) := nn >= 20 & cp_hi(kk, nn) < 0.025]; H0[, (paste0("rate_", v)) := fifelse(nn > 0, 100 * kk / nn, NA_real_)] }
H1[, `:=`(upow_orig = 100 * k / nd, upow_corr = 100 * (k - k_c) / nd)]
D <- merge(H0, H1[, .(dgp, scen, N, ntime, strategy, upow_orig, upow_corr)], by = c("dgp", "scen", "N", "ntime", "strategy")); D[, small := ntime <= 3 & N <= 25]
fwrite(D, "out/const2_dp.csv")
T5 <- D[, .(n_dp = .N, hi_orig = sum(hi_orig), hi_corr = sum(hi_corr), lo_orig = sum(lo_orig), lo_corr = sum(lo_corr), pu_orig = mean(pu_orig), pu_corr = mean(pu_corr),
            pu_small_orig = mean(pu_orig[small]), pu_small_corr = mean(pu_corr[small]), upow_orig = mean(upow_orig), upow_corr = mean(upow_corr),
            usable_const_H0 = sum(n_c)), by = strategy]
T5 <- merge(T5, H1[, .(usable_const_H1 = sum(n_c)), by = strategy], by = "strategy")
fwrite(T5, "out/const2_table5.csv"); print(T5, digits = 5)
