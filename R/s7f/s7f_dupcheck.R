source(path.expand("~/MMM/R/s4/s4_common.R"))
FL <- list.files(file.path(REPO, "results_v4/S7d/rows"), pattern = "\\.rds$", full.names = TRUE)
rd <- function(f) { x <- readRDS(f)$rows; x[cor == "exch", .(f = basename(f), dgp, scen, N, ntime, label, sim_id, gee_ok, est, se_kc)] }
F <- rbindlist(parallel::mclapply(FL, rd, mc.cores = 16))
print(F[, .N, by = scen])
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
d <- F[, .N, by = K][N > 1]
cat("duplicated keys", nrow(d), "\n"); print(head(d))
if (nrow(d)) { x <- merge(F, d[1:3, ..K], by = K); print(x) }
print(F[, .(n = .N, files = uniqueN(f)), by = .(scen)])
