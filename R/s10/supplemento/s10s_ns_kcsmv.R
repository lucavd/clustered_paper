
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2)
D <- "/home/user/MMM/results_v4/S10/supplemento"; K <- c("dgp","scen","N","ntime","label","sim_id")
R <- readRDS(file.path(D, "recs_s7d.rds"))[scen %in% c("M6","M7")]
R[, estd := generated %in% TRUE & !(design_degenerate %in% TRUE)][, const := estd & (sum_y == 0 | sum_y == n_obs)]
X <- readRDS("/home/user/MMM/results_v4/S6b/ns/rows_ns.rds")[method == "gee_kc"]
X <- merge(X, R[estd == TRUE, c(K, "const"), with = FALSE], by = K)
X[, use := converged %in% TRUE & is.finite(estimate) & is.finite(se) & is.finite(p) & abs(estimate) < 10 & se > 0 & se < 10 & !const]
X[, rt := abs(estimate / se) > qt(0.975, N - 4L)]
A <- X[label == "manual", .(n = sum(use), k = sum(use & rt)), by = .(dgp, scen, N, ntime)]
fwrite(A, "/home/user/MMM/results_v4/S10/supplemento/s10s_ns_kcsmv.csv"); cat(nrow(A), "\n")
