
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R"); src <- readLines("/home/user/MMM/R/s10/tabelle/s10tab_03_other.R")
eval(parse(text = src[3:grep("^run <- function", src) - 1]))
print(one_ns(file.path(REPO, "results_v4/S6b/ns/rows/GenBinaryY_m6_N10_T2_manual_c01.rds")))
print(one_cop(file.path(REPO, "results_v4/S7e/rows/copula_m3_N10_T10_manual_c01.rds"))$dp)
print(one_cop(file.path(REPO, "results_v4/S7e/rows/copula_m1_N105_T10_manual_c01.rds"))$dp[1:6])
z <- one_flip(file.path(REPO, "results_v4/S7d2/rows/copula_m1_N10_T10_manual_c01.rds")); print(z$dp); print(z$why)
