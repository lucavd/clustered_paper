# =============================================================================
# s6_figdata.R -- S6: data behind the nested-loop plot and the zip plot.
# Inputs : tables/V4_perf_primary.csv (S4), results_v4/S4/rows_v4.rds, tables/V4_S4_truth_by_cell.csv,
#          tables/V4_S2_degenerate_designs.csv
# Outputs: tables/V4_S6_nestedloop.csv  one row per design point x method: H0 rejection rate (own set,
#                                       usability U0) with its MCSE, main grid and M1c
#          tables/V4_S6_zip.csv         per-replicate estimate, SE and true value, H1 arm, GLMM-native M1,
#                                       T = 4, N in {10, 45, 145}, the seven Wald methods
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
P <- fread(file.path(TABDIR, "V4_perf_primary.csv"))
NL <- P[arm == "H0" & set == "own" & usab == "U0" & method %in% PRIM8,
        .(dgp, scen, N, ntime, method, n_obs, n_den, rej = 100 * rej, mcse_rej = 100 * mcse_rej)]
setorder(NL, dgp, scen, ntime, N, method)
fwrite(NL, file.path(TABDIR, "V4_S6_nestedloop.csv"))
cat(sprintf("nested loop: %d rows, %d design points, %d with n_obs = 0\n", nrow(NL), uniqueN(NL[, .(dgp, scen, N, ntime)]), NL[n_obs == 0, .N]))
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[!is.na(term) & dgp == "GLMM-native" & scen == "M1" & ntime == 4 & N %in% c(10, 45, 145) & label == "mmlb" & method %in% WALD7,
       .(dgp, scen, N, ntime, label, method, sim_id, estimate, se, p, U0)]
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[model == 1, .(N, ntime, sim_id)]
V[, deg := FALSE]; V[DEG, on = c("N", "ntime", "sim_id"), deg := TRUE]
TT <- fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))[dgp == "GLMM-native" & scen == "M1" & ntime == 4 & label == "mmlb"]
stopifnot(nrow(TT) == 1)
V[, truth := fifelse(method %in% MARGINAL_METHODS, TT$truth_marginal, TT$truth_conditional)]
V[, usable := U0 & !deg]
setorder(V, N, method, sim_id)
fwrite(V[, .(N, ntime, method, sim_id, estimate, se, truth, usable)], file.path(TABDIR, "V4_S6_zip.csv"))
print(V[, .(n = .N, usable = sum(usable), truth = truth[1]), by = .(N, method)])
