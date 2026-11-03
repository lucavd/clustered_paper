
source(path.expand("~/MMM/R/s7e/s7e_common.R")); data.table::setDTthreads(8L)
TK <- s7e_tasks()[label == "manual"]
L <- parallel::mclapply(TK$file, function(f) { r <- readRDS(f); r$study[method == "glmer" & term %in% "time:trt", .(scen, N, ntime, sim_id, estimate, se, converged)] }, mc.cores = 12L)
G <- rbindlist(L); rm(L)
RECS <- rbindlist(parallel::mclapply(TK$file, function(f) readRDS(f)$recs[, .(scen, N, ntime, sim_id, design_degenerate)], mc.cores = 12L))
G <- merge(G, RECS, by = c("scen", "N", "ntime", "sim_id"))[design_degenerate == FALSE & converged & is.finite(estimate) & abs(estimate) < 10 & is.finite(se) & se > 0 & se < 10]
S <- G[, .(n = .N, mean_est = mean(estimate), empse = sd(estimate), modse = sqrt(mean(se^2)), rob_ratio = median(se) / mad(estimate), rej_z = 100 * mean(abs(estimate / se) > qnorm(0.975))), by = .(scen, N, ntime)][, ratio := modse / empse]
S[, rate_implied := 100 * 2 * pnorm(-qnorm(0.975) * ratio)]
fwrite(S, file.path(TABDIR, "V4_S7e_glmm_se_ratio.csv"))
print(S[N >= 85, .(ratio = mean(ratio), rob_ratio = mean(rob_ratio), rej_z = mean(rej_z), rate_implied = mean(rate_implied), mean_est = mean(mean_est), bias_rel = mean(mean_est / empse)), by = .(scen, ntime)][order(scen, ntime)], digits = 3)
