
source(path.expand("~/MMM/R/s7e/s7e_common.R")); data.table::setDTthreads(8L)
TK <- s7e_tasks()[label == "manual" & scen %in% c("M3", "M4", "M5") & ntime %in% c(3L, 4L, 6L)]
G <- rbindlist(parallel::mclapply(TK$file, function(f) { r <- readRDS(f)
  merge(r$study[method == "glmer" & term %in% "time:trt", .(scen, N, ntime, sim_id, estimate, se, converged, status)], r$recs[, .(sim_id, design_degenerate)], by = "sim_id") }, mc.cores = 12L))
G <- G[design_degenerate == FALSE & converged & is.finite(estimate) & abs(estimate) < 10 & is.finite(se) & se > 0 & se < 10]
G[, singular := grepl("singular", status)]
X <- G[, .(n = .N, rej = 100 * mean(2 * pt(-abs(estimate / se), N - 4L - (scen == "M5")) < 0.05), med_se = median(se)), by = .(scen, ntime, singular)]
X[, share := 100 * n / sum(n), by = .(scen, ntime)]
print(X[order(scen, ntime, singular)], digits = 3); fwrite(X, file.path(TABDIR, "V4_S7e_glmm_singular.csv"))
