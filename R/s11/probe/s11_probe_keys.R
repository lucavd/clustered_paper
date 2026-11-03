# s11_probe_keys.R -- S11 read-only probe: datasets where the D-90 rule gives no result for the MMM with the check.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table)}); setDTthreads(4L); setwd("/home/user/MMM")
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
t0 <- Sys.time()
A <- readRDS("results_v4/S6cc/cc_all.rds"); A[, src := "main"]
fn <- list.files("results_v4/S6b/ns/rows", "\\.rds$", full.names = TRUE)
B <- rbindlist(parallel::mclapply(fn, function(f) readRDS(f)$cc, mc.cores = 8L), fill = TRUE); B[, src := "m67"]
fe <- list.files("results_v4/S7e/rows", "\\.rds$", full.names = TRUE)
C <- rbindlist(parallel::mclapply(fe, function(f) readRDS(f)$cc, mc.cores = 8L), fill = TRUE); C[, src := "copula"]
cols <- c(K6, "src", "A_ok", "A_has_term", "A_code", "A_slope_min", "bnd", "B_ok", "B_est", "B_se", "B_code", "used")
X <- rbind(A[, ..cols], B[, ..cols], C[, ..cols])
X[, sel := A_ok %in% TRUE & A_has_term %in% TRUE]
X[, bnd2 := sel & (!(A_code %in% 1:2) | (is.finite(A_slope_min) & A_slope_min < 1e-3))]
X[, Bconv := B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se) & B_code %in% 1:2]
X[, key := bnd2 & !Bconv]
print(X[, .(rows = .N, fitA = sum(sel), bnd_stored = sum(bnd %in% TRUE), bnd_recomp = sum(bnd2), bnd_mismatch = sum(sel & (bnd %in% TRUE) != bnd2),
            usedB = sum(used %in% "B"), usedB_notconv = sum(used %in% "B" & !Bconv), bnd_Bfail_keptA = sum(bnd2 & used %in% "A"), keys = sum(key)), by = src])
print(X[key == TRUE, .N, by = .(src, label)][order(src, label)])
cat("B codes of used refits, main grid:\n"); print(X[src == "main" & used %in% "B", .N, by = B_code][order(B_code)])
cat("duplicated keys:", anyDuplicated(X[, ..K6]), "\n")
cat("elapsed", format(Sys.time() - t0), "\n")
