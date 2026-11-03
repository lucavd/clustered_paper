# =============================================================================
# s11_cv_keys.R -- S11 (D-90, D-91): keys of the datasets where the MMM with the convergence check gives no result:
# the study fit (fit A) is at the boundary (nlm code not 1 or 2, or min exp(a_k) < 1e-3) and the fit without slope
# (fit B) does not converge (code not 1 or 2) or has no finite estimate and SE. Sources: results_v4/S6cc/cc_all.rds
# (main grid M3-M5), results_v4/S6b/ns/rows/*.rds $cc (M6-M7), results_v4/S7e/rows/*.rds $cc (Copula-AR1 M3-M5).
# Read-only on the sources. Output: results_v4/S11/cv_keys.rds (K6 + src + why), tables/V8_S11_cv_keys_summary.csv.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
A <- readRDS("results_v4/S6cc/cc_all.rds"); A[, src := "main"]
B <- rbindlist(parallel::mclapply(list.files("results_v4/S6b/ns/rows", "\\.rds$", full.names = TRUE), function(f) readRDS(f)$cc, mc.cores = 8L), fill = TRUE); B[, src := "m67"]
C <- rbindlist(parallel::mclapply(list.files("results_v4/S7e/rows", "\\.rds$", full.names = TRUE), function(f) readRDS(f)$cc, mc.cores = 8L), fill = TRUE); C[, src := "copula"]
cols <- c(K6, "src", "A_ok", "A_has_term", "A_code", "A_slope_min", "bnd", "B_ok", "B_est", "B_se", "B_code", "used")
X <- rbind(A[, ..cols], B[, ..cols], C[, ..cols]); rm(A, B, C)
stopifnot(!anyDuplicated(X, by = K6))
X[, fitA := A_ok %in% TRUE & A_has_term %in% TRUE]
X[, bnd2 := fitA & (!(A_code %in% 1:2) | (is.finite(A_slope_min) & A_slope_min < 1e-3))]
stopifnot(X[fitA == TRUE, all((bnd %in% TRUE) == bnd2)])                 # boundary as stored
X[, Bfin := B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se)]
X[, Bconv := Bfin & B_code %in% 1:2]
stopifnot(X[bnd2 & Bfin, all(used == "B")], X[bnd2 & !Bfin, all(used == "A")])  # stored strategy = variant (i)
KEYS <- X[bnd2 & !Bconv, c(K6, "src"), with = FALSE][, why := X[bnd2 & !Bconv, fifelse(Bfin, "refit not converged", "refit not finite")]]
setkeyv(KEYS, K6)
SUM <- X[, .(fits = sum(fitA), boundary = sum(bnd2), refit_used_v_i = sum(bnd2 & Bfin), refit_converged = sum(bnd2 & Bconv),
             no_result_not_converged = sum(bnd2 & Bfin & !Bconv), no_result_not_finite = sum(bnd2 & !Bfin), keys = sum(bnd2 & !Bconv)), by = .(src, label)]
print(SUM); fwrite(SUM, "tables/V8_S11_cv_keys_summary.csv")
dir.create("results_v4/S11", showWarnings = FALSE); saveRDS(KEYS, "results_v4/S11/cv_keys.rds")
cat("keys written:", nrow(KEYS), "\n")
