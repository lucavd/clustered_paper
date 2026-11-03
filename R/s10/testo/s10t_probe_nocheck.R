.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L); setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
cat("file:", file.exists("results_v4/S4/rows_v4.rds"), file.size("results_v4/S4/rows_v4.rds"), "\n")
S0 <- readRDS("results_v4/S4/rows_v4.rds")[method == "mmm" & scen %in% c("M1", "M2", "M3", "M4", "M5", "M1c")]
S0 <- S0[is.na(term) | term == "time:trt", c(K6, "term", "estimate", "se", "p", "converged", "status", "U0"), with = FALSE]
qs_save(S0, file.path(OUT, "x_study_nocheck_mmm.qs2"))
cat("original study mmm rows", nrow(S0), "\n")
C6 <- qs_read(file.path(OUT, "x_s6cc.qs2"))
Z <- merge(S0[scen %in% c("M3", "M4", "M5")], C6[, c(K6, "A_ok", "A_est", "A_se", "A_code", "bnd", "used", "cc_est"), with = FALSE], by = K6)
Z[, same := (abs(estimate - A_est) <= 1e-8) | (is.na(estimate) & is.na(A_est))]
cat("M3-M5 original study mmm vs S6cc fit A: pairs", nrow(Z), " same estimate", Z[same == TRUE, .N], " differ", Z[same == FALSE, .N], "\n")
print(Z[same == FALSE, .N, by = .(bnd, used, A_code)][order(-N)][1:10])
print(summary(Z[same == FALSE, abs(estimate - A_est)]))
V <- qs_read(file.path(OUT, "x_study.qs2"))[method == "mmm" & scen %in% c("M1", "M2", "M1c"), c(K6, "estimate"), with = FALSE]
W <- merge(S0[scen %in% c("M1", "M2", "M1c"), c(K6, "estimate"), with = FALSE], V, by = K6, suffixes = c("_orig", "_cc"))
cat("M1, M2, M1c: original vs variant_cc estimate differ:", W[!((abs(estimate_orig - estimate_cc) <= 1e-10) | (is.na(estimate_orig) & is.na(estimate_cc))), .N], "of", nrow(W), "\n")
