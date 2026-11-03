# s6_mmm_modes.R -- decomposition of the MMM excess rejection by failure mode (post hoc).
# Modes from the nlm code of the study fit A: "ok" = code 1 or 4 (converged, or stopped at the iteration
# limit with beta already at the maximum: see s6_mmm_iterlim_sum.R); "displaced" = code 5;
# "SE collapse" = code 2 or 3. For each scenario x T x N: share of each mode, its rejection, and its
# contribution to the excess = share x (rejection - 5). For code 5, the displacement |A - B| in absolute
# terms and relative to the SE, by N. Input results_v4/S6/mmm_boundary_all.rds.
# Output: tables/V4_S6_mmm_modes.csv, V4_S6_mmm_modes_disp.csv
source(path.expand("~/MMM/R/s6/s6_common.R"))
X <- readRDS(file.path(REPO, "results_v4/S6/mmm_boundary_all.rds"))[useA == TRUE & model %in% c(1L, 3:5)]
X[, mode := fcase(A_code %in% c(1L, 4L), "ok", A_code == 5L, "displaced", A_code %in% c(2L, 3L), "se_collapse", default = "other")]
pr <- function(x) 100 * mean(x, na.rm = TRUE)
M <- X[, .(n = .N), by = .(model, ntime, N)]
D <- X[, .(k = .N, rej = pr(rA)), by = .(model, ntime, N, mode)]
D <- merge(D, M, by = c("model", "ntime", "N")); D[, share := 100 * k / n][, contrib := share / 100 * (rej - 5)]
fwrite(D[order(model, ntime, N, mode)], file.path(TABDIR, "V4_S6_mmm_modes.csv"))
Q <- X[mode == "displaced" & useB == TRUE & model %in% 3:5, .(n = .N, med_abs_disp = median(abs(A_est - B_est)), q90_abs_disp = quantile(abs(A_est - B_est), .9),
        med_seB = median(B_se_mod), q90_disp_over_se = quantile(abs(A_est - B_est) / B_se_mod, .9), rejA = pr(rA), rejB = pr(rB)), by = .(model, N)]
fwrite(Q[order(model, N)], file.path(TABDIR, "V4_S6_mmm_modes_disp.csv"))
options(width = 220)
T8 <- D[model %in% 3:5 & N >= 25, .(share = weighted.mean(share, n), rej = weighted.mean(rej, k), contrib = weighted.mean(contrib, n)), by = .(model, ntime, mode)]
print(dcast(T8, model + ntime ~ mode, value.var = c("share", "rej", "contrib")), digits = 2)
NN <- D[model %in% 3:5 & ntime >= 8, .(share = weighted.mean(share, n), rej = weighted.mean(rej, k), contrib = weighted.mean(contrib, n)), by = .(model, N, mode)]
print(dcast(NN, model + N ~ mode, value.var = c("share", "rej", "contrib")), digits = 2)
print(Q[order(model, N)], digits = 3)
