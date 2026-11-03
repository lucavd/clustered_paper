.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L); setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
C6 <- qs_read(file.path(OUT, "x_s6cc.qs2")); C7 <- qs_read(file.path(OUT, "x_nscc.qs2")); C8 <- qs_read(file.path(OUT, "x_cop_cc.qs2"))
cols <- c(K6, "A_ok", "A_est", "A_se", "A_p", "A_code", "A_slope_min", "bnd", "B_ok", "B_est", "B_se", "B_p", "B_code", "used", "cc_est", "cc_se")
CC <- rbind(C6[, ..cols][, src := "main"], C7[, ..cols][, src := "m67"], C8[, ..cols][, src := "copula"])
X <- CC[bnd %in% TRUE & B_ok %in% TRUE & used == "A"]
cat("bnd & B_ok & used A:", nrow(X), "\n"); print(X[, .N, by = .(src, B_code, Bfin = is.finite(B_est) & is.finite(B_se))])
print(X[, .(A_code = paste(names(table(A_code)), table(A_code), collapse = " ")), by = src])
print(head(X[, .(src, dgp, scen, N, ntime, A_code, A_slope_min, A_est, B_est, B_se, B_code, cc_est, cc_se)], 8))
Y <- CC[bnd %in% TRUE & B_ok %in% TRUE & used == "B"]; print(Y[, .N, by = .(src, B_code)])
cat("cc_est == A_est in used A:", X[, sum(abs(cc_est - A_est) < 1e-12, na.rm = TRUE)], "\n")
S <- qs_read(file.path(OUT, "x_study.qs2"))[method == "mmm" & scen %in% c("M3", "M4", "M5", "M6", "M7"), c(K6, "estimate", "se"), with = FALSE]
Z <- merge(S, CC[src != "copula"], by = K6)
Z[, d_cc := abs(estimate - cc_est)]; Z[, d_A := abs(estimate - A_est)]; Z[, d_B := abs(estimate - B_est)]
cat("stored mmm vs stored cc_est: differ", Z[!(d_cc <= 1e-10 | (is.na(estimate) & is.na(cc_est))), .N], "\n")
W <- Z[!(d_cc <= 1e-10 | (is.na(estimate) & is.na(cc_est)))]
print(W[, .N, by = .(src, used, est_na = is.na(estimate), cc_na = is.na(cc_est), se_na = is.na(se))])
print(summary(W$d_cc)); print(head(W[, .(scen, N, ntime, used, estimate, se, cc_est, cc_se, A_est, B_est)], 6))
