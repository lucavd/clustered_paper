# =============================================================================
# s11_cc_cells.R -- S11 (D-90, D-91): tables/V4_S6cc_cells.csv (check vs study fit, per design point and arm, read by
# R/s7f/s7f_cpflags.R for Section 3.2) recomputed with the MMM with the check taken from the variant cv instead of cc.
# Same computation as R/s6/s6cc_compare.R lines 11-16 and 48-58 (the S6 predictions and integrity checks of that
# script are not rerun). Gate: with the variant cc the table stored at the close of S10 (commit 7a68326) is reproduced exactly.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
stopifnot(!nzchar(VARIANT))
VS <- readRDS(file.path(REPO, "results_v4/S4/rows_v4.rds"))
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
prep <- function(V) {
  X <- V[method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term), .(dgp, scen, N, ntime, label, sim_id, p, U0)]
  X[, model := as.integer(substr(scen, 2, 2))]; X[, deg := FALSE]; X[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
  X[, u := U0 & !deg]; X[, r := u & fifelse(is.na(p), FALSE, p < 0.05)]; X[, .(dgp, scen, model, N, ntime, label, sim_id, u, r)] }
S <- prep(VS); rm(VS); invisible(gc(verbose = FALSE))
K <- c("dgp", "scen", "model", "N", "ntime", "label", "sim_id")
cells <- function(C) {
  W <- merge(S, C, by = K, suffixes = c("_s", "_c")); stopifnot(nrow(W) == nrow(S), nrow(W) == nrow(C))
  H0 <- W[label == "manual", .(n_c = sum(u_c), k_c = sum(r_c), n_s = sum(u_s), k_s = sum(r_s)), by = .(dgp, scen, N, ntime)]
  H0[, `:=`(rate_c = 100 * k_c / n_c, rate_s = 100 * k_s / n_s)]
  H0[, `:=`(hi_c = rate_c - 200 * sqrt(rate_c / 100 * (1 - rate_c / 100) / n_c) > 7.5, hi_s = rate_s - 200 * sqrt(rate_s / 100 * (1 - rate_s / 100) / n_s) > 7.5)]
  H0[, `:=`(lo_c = rate_c + 200 * sqrt(rate_c / 100 * (1 - rate_c / 100) / n_c) < 2.5, lo_s = rate_s + 200 * sqrt(rate_s / 100 * (1 - rate_s / 100) / n_s) < 2.5)]
  H1 <- W[label == "mmlb", .(n_c = sum(u_c), k_c = sum(r_c), n_s = sum(u_s), k_s = sum(r_s)), by = .(dgp, scen, N, ntime)]
  rbind(H0[, .(arm = "H0", dgp, scen, N, ntime, n_c, k_c, n_s, k_s, hi_c, hi_s, lo_c, lo_s)], H1[, .(arm = "H1", dgp, scen, N, ntime, n_c, k_c, n_s, k_s)], fill = TRUE)
}
OLD <- fread(cmd = "git -C /home/user/MMM show 7a68326:tables/V4_S6cc_cells.csv")   # stored table at the close of S10 (before S11)
VC <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds")); CEcc <- cells(prep(VC)); rm(VC); invisible(gc(verbose = FALSE))
fwrite(CEcc, "/tmp/s11_cells_cc.csv"); R1 <- fread("/tmp/s11_cells_cc.csv")
rep_ok <- isTRUE(all.equal(R1, OLD, check.attributes = FALSE))
cat("reproduction of the stored table with the variant cc:", rep_ok, "\n"); stopifnot(rep_ok)
VV <- readRDS(file.path(REPO, "results_v4/variant_cv/S4/rows_v4.rds")); CE <- cells(prep(VV)); rm(VV); invisible(gc(verbose = FALSE))
same_s <- CE[, .(arm, dgp, scen, N, ntime, n_s, k_s)]; same_s0 <- CEcc[, .(arm, dgp, scen, N, ntime, n_s, k_s)]
stopifnot(isTRUE(all.equal(same_s, same_s0)))   # the study side does not change
cat("H0 points clearly above 7.5%: check (cc)", CEcc[arm == "H0", sum(hi_c, na.rm = TRUE)], "check (cv)", CE[arm == "H0", sum(hi_c, na.rm = TRUE)],
    "study", CE[arm == "H0", sum(hi_s, na.rm = TRUE)], "\n")
fwrite(CE, file.path(TABDIR, "V4_S6cc_cells.csv"))
fwrite(data.table(check = c("cc reproduces stored cells", "study side unchanged", "hi_c cc", "hi_c cv", "lo_c cc", "lo_c cv"),
                  value = c(rep_ok, TRUE, CEcc[arm == "H0", sum(hi_c, na.rm = TRUE)], CE[arm == "H0", sum(hi_c, na.rm = TRUE)],
                            CEcc[arm == "H0", sum(lo_c, na.rm = TRUE)], CE[arm == "H0", sum(lo_c, na.rm = TRUE)])), file.path(TABDIR, "V8_S11_cc_cells_gate.csv"))
