
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
x <- as.data.table(readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))); o <- own_cc(x)
d <- cbind(x, o[, .(used_own, bnd_own)])
a <- d[used_own == "B" & used == "A"]
print(a[, .N, by = .(A_code, slope_lt = A_slope_min < 0.001, B_ok, bnd, A_has_term)][order(-N)])
print(summary(a$A_slope_min)); print(head(a[, .(A_code, A_slope_min, A_alpha, A_alpha_names, bnd, B_ok, B_code)], 5))
b <- d[is.na(used)]; print(b[, .N, by = .(A_ok, A_has_term, A_code_na = is.na(A_code), B_ok, bnd, scen)])
print(head(b[A_ok == TRUE, .(A_est, A_se, A_p, A_code, A_status, A_slope_min, cc_est, cc_status)], 5))
print(x[, .N, by = .(bnd, bnd_own = o$bnd_own)])
