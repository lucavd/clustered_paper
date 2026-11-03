# =============================================================================
# s6_mmm_diag.R -- S6 follow-up of the PI's question: why the marginalized model's departure
# from 5% rises with series length (Figure 1b) and its share of clearly anticonservative design
# points rises with N (Figure 1c). Null arm, own set, usability U0, from tables/V4_perf_primary.csv:
# per design point the rejection rate, bias of the estimate (true value 0 under H0, marginal),
# empirical SE, model SE and their ratio, and the rejection rate a Wald test would have if the
# estimate were normal with that bias and empirical SE and were divided by the model SE:
#   P(|Z| > 1.96), Z ~ N(bias / modSE, (empSE / modSE)^2)
# split into "bias only" (empSE = modSE) and "SE only" (bias = 0).
# Output: tables/V4_S6_mmm_diag.csv (design points), tables/V4_S6_mmm_diag_summary.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
P <- fread(file.path(TABDIR, "V4_perf_primary.csv"))[arm == "H0" & set == "own" & usab == "U0" & method %in% c("mmm", "gee_fg", "gee_md")]
pr <- function(mu, s) pnorm(-1.959964, mu, s) + pnorm(1.959964, mu, s, lower.tail = FALSE)
P[, `:=`(rate = 100 * rej, z_bias = bias / modse, se_ratio = modse / empse)]
P[, `:=`(pred = 100 * pr(z_bias, 1 / se_ratio), pred_bias_only = 100 * pr(z_bias, 1), pred_se_only = 100 * pr(0, 1 / se_ratio))]
X <- P[, .(dgp, scen, N, ntime, method, n_obs, rate, mcse_rej = 100 * mcse_rej, bias, mcse_bias, empse, modse, se_ratio, z_bias, pred, pred_bias_only, pred_se_only)]
fwrite(X, file.path(TABDIR, "V4_S6_mmm_diag.csv"))
S <- X[n_obs >= 20, .(n_dp = .N, rate = mean(rate), dev = mean(abs(rate - 5)), z_bias = mean(z_bias), abs_z_bias = mean(abs(z_bias)), se_ratio = mean(se_ratio),
                      pred = mean(pred), pred_bias_only = mean(pred_bias_only), pred_se_only = mean(pred_se_only)), by = .(method, dgp, scen, ntime)]
setorder(S, method, scen, dgp, ntime); fwrite(S, file.path(TABDIR, "V4_S6_mmm_diag_summary.csv"))
print(dcast(S[method == "mmm"], scen + dgp ~ ntime, value.var = "dev"), digits = 2)
print(dcast(S[method == "mmm"], scen + dgp ~ ntime, value.var = "z_bias"), digits = 2)
print(dcast(S[method == "mmm"], scen + dgp ~ ntime, value.var = "se_ratio"), digits = 3)
