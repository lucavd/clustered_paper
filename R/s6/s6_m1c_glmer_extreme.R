# =============================================================================
# s6_m1c_glmer_extreme.R -- S6 follow-up: the size-adjusted critical p-value of glmer in M1c
# is ~0 in some cells (GLMM-native, T = 6), which drives its size-adjusted power to 0-2%.
# Counts usable null-arm replicates with an extreme Wald p-value (p < 1e-10), per cell and method,
# and describes them (estimate, SE, status). Output: tables/V4_S6_m1c_glmer_extreme.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[!is.na(term) & method %in% c("glmer", "glmer_bobyqa", "shrink", "mmm", "gee_fg") & U0 == TRUE & design_degenerate %in% FALSE]
V[, grid := fifelse(scen == "M1c", "M1c", "main")]
X <- V[, .(n = .N, n_extreme = sum(p < 1e-10), pct_extreme = 100 * mean(p < 1e-10),
           med_abs_est_extreme = median(abs(estimate[p < 1e-10])), med_se_extreme = median(se[p < 1e-10]),
           med_se_all = median(se), warn_extreme = 100 * mean(nzchar(status[p < 1e-10]) & !is.na(status[p < 1e-10]))),
       by = .(grid, dgp, label, method, ntime)]
setorder(X, grid, dgp, method, label, ntime); fwrite(X, file.path(TABDIR, "V4_S6_m1c_glmer_extreme.csv"))
print(X[method %in% c("glmer", "glmer_bobyqa") & (grid == "M1c" | ntime >= 6)], nrows = 200)
ex <- V[grid == "M1c" & method == "glmer" & p < 1e-10][order(se)]
print(head(ex[, .(dgp, N, ntime, label, sim_id, estimate, se, p, status = substr(status, 1, 80))], 15))
print(ex[, .N, by = .(status = substr(status, 1, 60))][order(-N)][1:10])
