# =============================================================================
# s2_f_pilot_analysis.R -- S2(f), Q6.3/Q6.4: design identity and paired
# comparison corrected H0 (results_v4/pilot_S2) vs original H0 (results_full).
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
PILOT <- file.path(REPO, "results_v4", "pilot_S2")
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink"); LRT <- "glmer_lrt"
fs <- list.files(PILOT, pattern = "manualv4_c[0-9]+[.]rds$", full.names = TRUE)
cat("pilot chunks:", length(fs), "\n"); stopifnot(length(fs) == 120L)
kre <- "^GenBinaryY_m([3-5])_N([0-9]+)_T([0-9]+)_manualv4_c([0-9]+)[.]rds$"
rd <- function(f, keepdes = FALSE) {
  x <- readRDS(f); b <- basename(f)
  des <- as.data.table(attr(x, "design_check"))
  x <- as.data.table(x)[is.na(term) | term == "time:trt", .(method, term, estimate, se, p, converged, status, sim_id = as.integer(sim_id))]
  x[, `:=`(model = as.integer(sub(kre, "\\1", b)), N = as.integer(sub(kre, "\\2", b)), ntime = as.integer(sub(kre, "\\3", b)))]
  des[, `:=`(model = x$model[1], N = x$N[1], ntime = x$ntime[1])]
  list(x = x, des = des)
}
L <- lapply(fs, rd)
NEW <- rbindlist(lapply(L, `[[`, "x")); DES <- rbindlist(lapply(L, `[[`, "des"))
cat(sprintf("Q6.3 design identical (id, time, trt, dose) corrected vs original H0: %d/%d datasets\n", sum(DES$design_identical), nrow(DES)))
cat("Q6.3 share of identical binY values, by cell (informative):\n"); print(DES[, .(y_same = mean(y_identical_share)), by = .(model, N, ntime)])
A <- readRDS(file.path(OUTDIR, "rows_timetrt.rds"))[dgp == "GenBinaryY" & label == "manual" & model %in% 3:5]
cells <- unique(NEW[, .(model, N, ntime)])
OLD <- A[cells, on = .(model, N, ntime), nomatch = 0][, .(method, term, estimate, se, p, converged, status, sim_id, model, N, ntime)]
use <- function(D) {
  D <- D[method != "gee_ss"]
  D[, U := fifelse(method %in% WALD, converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p),
                   converged & is.finite(p))]
  D[is.na(U), U := FALSE]; D[, R := U & p < 0.05]; D[, .(model, N, ntime, method, sim_id, U, R)]
}
O <- use(OLD); Nw <- use(NEW)
grid <- CJ(model = 3:5, N_ntime = 1:4, method = c(WALD, LRT), sim_id = 1:1000)[, `:=`(N = c(10, 25, 65, 145)[N_ntime], ntime = c(2, 4, 6, 10)[N_ntime])][, N_ntime := NULL]
M <- merge(grid, O, by = c("model", "N", "ntime", "method", "sim_id"), all.x = TRUE)
M <- merge(M, Nw, by = c("model", "N", "ntime", "method", "sim_id"), all.x = TRUE, suffixes = c("_old", "_new"))
for (v in c("U_old", "U_new", "R_old", "R_new")) M[is.na(get(v)), (v) := FALSE]
S <- M[, {
  du <- as.numeric(U_new) - as.numeric(U_old); b <- U_new & U_old; dr <- as.numeric(R_new[b]) - as.numeric(R_old[b])
  .(n = .N, usable_old = mean(U_old), usable_new = mean(U_new), d_usable = mean(du), mcse_d_usable = sd(du) / sqrt(.N),
    t1_old = sum(R_old) / max(sum(U_old), 1), t1_new = sum(R_new) / max(sum(U_new), 1),
    mcse_t1_new = mcse_p(sum(R_new) / max(sum(U_new), 1), max(sum(U_new), 1)),
    n_both = sum(b), d_t1_paired = if (sum(b)) mean(dr) else NA_real_, mcse_d_t1 = if (sum(b) > 1) sd(dr) / sqrt(sum(b)) else NA_real_)
}, by = .(model, N, ntime, method)]
z3 <- function(d, s) fifelse(is.na(d), NA, fifelse(s == 0, d == 0, abs(d) <= 3 * s))
S[, `:=`(pass_usable = z3(d_usable, mcse_d_usable), pass_t1 = z3(d_t1_paired, mcse_d_t1))]
fwrite(S, file.path(TABDIR, "V4_S2_pilot_h0.csv"))
options(width = 220); print(S, digits = 3)
k <- S[model %in% 3:4]
cat(sprintf("Q6.4 M3-M4: usable within 3 MCSE %d/%d; type I (paired) within 3 MCSE %d/%d (NA %d)\n",
            sum(k$pass_usable, na.rm = TRUE), sum(!is.na(k$pass_usable)), sum(k$pass_t1, na.rm = TRUE), sum(!is.na(k$pass_t1)), sum(is.na(k$pass_t1))))
cat("Q6.4 M3-M4 exceedances:\n"); print(k[pass_usable == FALSE | pass_t1 == FALSE])
cat("Q6.4 M5 (exploratory) exceedances:\n"); print(S[model == 5 & (pass_usable == FALSE | pass_t1 == FALSE)])
