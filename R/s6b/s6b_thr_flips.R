# s6b_thr_flips.R -- S6b part A, explanation of PA1-PA3: in the fits that change with the threshold, how much do
# the estimate and the test decision move between the study fit (A) and the fit without slope (B)?
# Output: tables/V4_S6b_thr_flips.csv
source(path.expand("~/MMM/R/s4/s4_common.R"))
CC <- readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))[A_ok %in% TRUE & A_has_term %in% TRUE]
TH <- rbindlist(lapply(list.files(file.path(REPO, "results_v4/S6b/thr"), "^rows_.*\\.rds$", full.names = TRUE), readRDS), fill = TRUE)[set == "a"]
conv <- CC$A_code %in% c(1L, 2L)
s4 <- CC[conv & A_slope_min >= 1e-4 & A_slope_min < 1e-3 & used == "B", .(label, eA = A_est, sA = A_se, pA = A_p, eB = B_est, sB = B_se, pB = B_p)][, change := "1e-3 -> 1e-4: B back to A"]
s2 <- TH[B2_ok %in% TRUE & is.finite(B2_est) & is.finite(B2_se), .(label, eA = A_est, sA = A_se, pA = A_p, eB = B2_est, sB = B2_se, pB = B2_p)][, change := "1e-3 -> 1e-2: A to B"]
X <- rbind(s4, s2)
F <- X[, .(n = .N, med_abs_dz = median(abs(eA - eB) / sA, na.rm = TRUE), q99_abs_dz = quantile(abs(eA - eB) / sA, 0.99, na.rm = TRUE), med_se_ratio = median(sB / sA, na.rm = TRUE),
           n_na = sum(!is.finite(pA) | !is.finite(pB)), rej_A = sum(pA < 0.05, na.rm = TRUE), rej_B = sum(pB < 0.05, na.rm = TRUE), flip_A_only = sum(pA < 0.05 & !(pB < 0.05), na.rm = TRUE), flip_B_only = sum(pB < 0.05 & !(pA < 0.05), na.rm = TRUE)), by = .(change, label)]
fwrite(F, file.path(TABDIR, "V4_S6b_thr_flips.csv")); print(F, digits = 3)
