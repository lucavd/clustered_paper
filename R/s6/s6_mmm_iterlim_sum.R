# s6_mmm_iterlim_sum.R -- result of the intervention test (s6_mmm_iterlim.R).
# For the non-converged fits of M3-M5, T >= 8, N = 55, 145: nlm code after 2000 iterations,
# rejection and SE ratio before (study fit, 100 iterations) and after; for the whole cell, the study
# fit with only the non-converged fits replaced. Output: tables/V4_S6_mmm_iterlim.csv, V4_S6_mmm_iterlim_codes.csv
source(path.expand("~/MMM/R/s6/s6_common.R"))
X <- readRDS(file.path(REPO, "results_v4/S6/mmm_boundary_all.rds"))[useA == TRUE & model %in% 3:5 & ntime >= 8 & N %in% c(55L, 145L)]
L <- rbindlist(lapply(c("GenBinaryY", "GLMM-native", "GLMM-AR1"), function(g) readRDS(file.path(REPO, sprintf("results_v4/S6/mmm_iterlim_%s.rds", g)))), fill = TRUE)
X <- merge(X, L, by = c("dgp", "model", "ntime", "N", "i"), all.x = TRUE)
z <- qnorm(0.975); pr <- function(x) 100 * mean(x, na.rm = TRUE)
X[, L_use := L_ok %in% TRUE & is.finite(L_est) & is.finite(L_se_mod)]
X[, L_conv := L_use & L_code %in% c(1L, 2L)]
X[, rL := abs(L_est / L_se_mod) > z]
X[, fix_est := fifelse(nonconv & L_use, L_est, A_est)][, fix_se := fifelse(nonconv & L_use, L_se_mod, A_se_mod)][, rFix := abs(fix_est / fix_se) > z]
S1 <- X[nonconv == TRUE, .(n = .N, n_refit = sum(L_use), pct_conv_after = pr(L_conv), med_iter_after = as.numeric(median(L_niter, na.rm = TRUE)),
        rej_before = pr(rA), rej_after = pr(rL[L_use]), ratio_before = sqrt(mean(A_se_mod^2)) / sd(A_est),
        ratio_after = sqrt(mean(L_se_mod[L_use]^2)) / sd(L_est[L_use]), rej_after_conv = pr(rL[L_conv]), n_conv_after = sum(L_conv),
        rej_B = pr(rB[useB])), by = .(model, N)]
S2 <- X[, .(n = .N, rej_study = pr(rA), rej_fixed = pr(rFix), rej_B = pr(rB[useB]), ratio_fixed = sqrt(mean(fix_se^2)) / sd(fix_est)), by = .(model, N)]
S <- merge(S1, S2, by = c("model", "N"), suffixes = c("_nonconv", "_cell"))
S[, mcse_fixed := 100 * sqrt(rej_fixed / 100 * (1 - rej_fixed / 100) / n_cell)]
fwrite(S[order(model, N)], file.path(TABDIR, "V4_S6_mmm_iterlim.csv"))
CO <- X[nonconv == TRUE & L_use == TRUE, .N, by = .(model, code_before = A_code, code_after = L_code)][order(model, code_before, code_after)]
fwrite(CO, file.path(TABDIR, "V4_S6_mmm_iterlim_codes.csv"))
options(width = 250); print(S[order(model, N)], digits = 3); print(CO)
