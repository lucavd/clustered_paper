
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2)
D <- "/home/user/MMM/results_v4/S10/supplemento"; K <- c("dgp","scen","N","ntime","label","sim_id")
E <- readRDS(file.path(D, "study_s7e.rds"))[method == "glmer" & term == "time:trt" & scen %in% c("M1","M2")]
R <- readRDS(file.path(D, "recs_s7e.rds"))[scen %in% c("M1","M2")]
R[, estd := generated %in% TRUE & !(design_degenerate %in% TRUE)][, const := estd & (sum_y == 0 | sum_y == n_obs)]
X <- merge(E, R[estd == TRUE, c(K, "const"), with = FALSE], by = K)
X <- X[converged %in% TRUE & is.finite(estimate) & abs(estimate) < 10 & is.finite(se) & se > 0 & se < 10 & !const & N >= 85]
P <- X[, .(n = .N, sde = sd(estimate), rms = sqrt(mean(se^2)), mse = mean(se)), by = .(scen, N, ntime, label)]
P[, `:=`(r_rms = rms / sde, r_mean = mse / sde)]
for (lab in c("manual", "mmlb")) for (v in c("r_rms", "r_mean")) {
  G <- P[label == lab, .(r = mean(get(v))), by = .(scen, ntime)]
  cat(sprintf("%-6s %-6s by(scen,T) mean over N: min %.4f max %.4f | per dp: min %.4f max %.4f\n", lab, v, min(G$r), max(G$r), min(P[label == lab][[v]]), max(P[label == lab][[v]])))
}
G <- P[label == "manual", .(r = mean(r_rms)), by = .(scen, ntime)][order(r)]; print(G)
cat("const among glmer-usable N>=85:", X[, sum(const)], "\n")
