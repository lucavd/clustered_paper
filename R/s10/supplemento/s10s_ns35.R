
# s10s_ns35.R -- S10 "supplemento": Table S35 SE ratio by intervention of the check, with the mean model SE.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(parallel)}); setDTthreads(2)
RP <- "/home/user/MMM/results_v4"; D <- file.path(RP, "S10/supplemento"); K <- c("dgp","scen","N","ntime","label","sim_id"); CELL <- K[1:5]
R <- readRDS(file.path(D, "recs_s7d.rds"))[scen %in% c("M6","M7")]
R[, estd := generated %in% TRUE & !(design_degenerate %in% TRUE)][, const := estd & (sum_y == 0 | sum_y == n_obs)]
f <- list.files(file.path(RP, "S6b/ns/rows"), full.names = TRUE)
X <- rbindlist(mclapply(f, function(x) readRDS(x)$cc, mc.cores = 2), fill = TRUE)
X <- merge(X, R[estd == TRUE, c(K, "const"), with = FALSE], by = K)
X[, use := is.finite(cc_p) & is.finite(cc_est) & is.finite(cc_se) & abs(cc_est) < 10 & cc_se > 0 & cc_se < 10 & !const][, w := used == "B"]
P <- X[use == TRUE, .(n = .N, n_w = sum(w), mse_w = mean(cc_se[w]), mse_wo = mean(cc_se[!w]), sde_w = sd(cc_est[w]), sde_wo = sd(cc_est[!w]), sde = sd(cc_est)), by = CELL]
fwrite(P, file.path(D, "s10s_ns35_sub.csv")); cat("NS35 OK\n")
