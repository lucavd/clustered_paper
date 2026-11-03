# =============================================================================
# s11_cc_shares.R -- S11 (D-90): how often the convergence check of the MMM intervened, per set of scenarios, arm
# and band of occasions: share of study fits declared at the boundary, replaced by the converged fit without slope,
# and left without a result (refit not converged or not finite). Denominator: study fits that returned the term.
# Reads the stored check fits (read-only). Output: tables/V8_S11_cc_shares.csv.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
A <- readRDS("results_v4/S6cc/cc_all.rds"); A[, set := "main M3-M5"]
B <- rbindlist(parallel::mclapply(list.files("results_v4/S6b/ns/rows", "\\.rds$", full.names = TRUE), function(f) readRDS(f)$cc, mc.cores = 8L), fill = TRUE); B[, set := "M6-M7"]
C <- rbindlist(parallel::mclapply(list.files("results_v4/S7e/rows", "\\.rds$", full.names = TRUE), function(f) readRDS(f)$cc, mc.cores = 8L), fill = TRUE); C[, set := "Copula-AR1 M3-M5"]
cols <- c("dgp", "scen", "N", "ntime", "label", "sim_id", "set", "A_ok", "A_has_term", "A_code", "A_slope_min", "B_ok", "B_est", "B_se", "B_code")
X <- rbind(A[, ..cols], B[, ..cols], C[, ..cols])[A_ok %in% TRUE & A_has_term %in% TRUE]
X[, `:=`(bnd = !(A_code %in% 1:2) | (is.finite(A_slope_min) & A_slope_min < 1e-3), nonconvA = !(A_code %in% 1:2),
         Bconv = B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se) & B_code %in% 1:2)]
X[, arm := fifelse(label == "manual", "H0", "H1")]; X[, band := fifelse(ntime <= 3, "T2-3", fifelse(ntime <= 6, "T4-6", "T8-10"))]
f <- function(D) D[, .(fits = .N, boundary = 100 * mean(bnd), nonconv_A = 100 * mean(nonconvA), replaced = 100 * mean(bnd & Bconv), no_result = 100 * mean(bnd & !Bconv),
                       no_result_of_boundary = 100 * sum(bnd & !Bconv) / sum(bnd))]
OUT <- rbind(X[, f(.SD), by = .(set, arm)][, `:=`(scen = "all", band = "all", dgp = "all")],
             X[, f(.SD), by = .(set, arm, band)][, `:=`(scen = "all", dgp = "all")],
             X[, f(.SD), by = .(set, arm, dgp, scen)][, band := "all"],
             X[, f(.SD), by = .(set, arm, dgp, scen, band)], use.names = TRUE)
setcolorder(OUT, c("set", "arm", "dgp", "scen", "band"))
fwrite(OUT, "tables/V8_S11_cc_shares.csv"); print(OUT[scen == "all"], digits = 3)
