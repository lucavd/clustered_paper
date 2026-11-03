# s6_mmm_mechanism.R -- why the MMM excess grows with T and N (post hoc, on results_v4/S6/mmm_boundary_all.rds).
# (1) rejection of the study fit (A) and of the fit without latent-SD slope (B) by scenario x T x N,
#     pooled over the three mechanisms; nlm code shares.
# (2) for each fit, the gap between A and B estimates in units of A's SE: if A stops before the maximum,
#     the gap should not shrink with N while the SE does, so gap/SE grows with N.
# (3) M1 cells outside the 2.5-7.5 band (P7 detail).
# Outputs: tables/V4_S6_mmm_mech_TN.csv, V4_S6_mmm_mech_gap.csv, V4_S6_mmm_mech_codes.csv, V4_S6_mmm_mech_m1.csv
source(path.expand("~/MMM/R/s6/s6_common.R"))
X <- readRDS(file.path(REPO, "results_v4/S6/mmm_boundary_all.rds"))[useA == TRUE]
pr <- function(x) 100 * mean(x, na.rm = TRUE)
TN <- X[, .(n = .N, rej_A = pr(rA), rej_B = pr(rB), rej_rule = pr(rRule), nonconv = pr(nonconv), code3 = pr(A_code == 3), code4 = pr(A_code == 4), code5 = pr(A_code == 5),
            se_ratio_A = sqrt(mean(A_se_mod^2)) / sd(A_est), se_ratio_B = sqrt(mean(B_se_mod[useB]^2)) / sd(B_est[useB])), by = .(model, ntime, N)]
TN[, mcse := 100 * sqrt(rej_A / 100 * (1 - rej_A / 100) / n)]
fwrite(TN[order(model, ntime, N)], file.path(TABDIR, "V4_S6_mmm_mech_TN.csv"))
G <- X[model %in% 3:5 & useB == TRUE, .(n = .N, med_gap = median(abs(A_est - B_est)), med_seA = median(A_se_mod), med_gap_over_se = median(abs(A_est - B_est) / A_se_mod),
        q90_gap_over_se = quantile(abs(A_est - B_est) / A_se_mod, 0.9)), by = .(model, conv = !nonconv, Tb, N)]
fwrite(G[order(model, conv, Tb, N)], file.path(TABDIR, "V4_S6_mmm_mech_gap.csv"))
CO <- X[, .N, by = .(model, Tb, code = A_code)][, pct := 100 * N / sum(N), by = .(model, Tb)]
fwrite(CO[order(model, Tb, code)], file.path(TABDIR, "V4_S6_mmm_mech_codes.csv"))
C <- fread(file.path(TABDIR, "V4_S6_mmm_bnd_cells.csv"))
M1 <- C[model == 1 & n_use >= 100 & (rej_all < 2.5 | rej_all > 7.5), .(dgp, ntime, N, n_use, rej_all, mcse_all, pct_nonconv)]
fwrite(M1, file.path(TABDIR, "V4_S6_mmm_mech_m1.csv"))
options(width = 250)
print(dcast(TN[model %in% 3:5], model + ntime ~ N, value.var = "rej_A"), digits = 3)
print(dcast(TN[model %in% 3:5], model + ntime ~ N, value.var = "rej_B"), digits = 3)
print(dcast(TN[model %in% 3:5], model + ntime ~ N, value.var = "nonconv"), digits = 3)
print(dcast(TN[model %in% 3:5], model + ntime ~ N, value.var = "se_ratio_A"), digits = 3)
print(G[Tb == "T>=8"], digits = 3)
print(CO[Tb == "T>=8"], digits = 3)
print(M1, digits = 3)
