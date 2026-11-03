# =============================================================================
# s6_dp_rates.R -- S6, rule evaluation: per design point and method, null arm, own set,
# usability U0: usable replicates n, non-degenerate replicates nd, rejections under the
# normal (k) and the t reference with N - p degrees of freedom (kt), from row b = 0 of the
# S5 aggregates (all 1,000 replicates once). Main grid and M1c.
# Output: tables/V4_S6_dp_rates.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
AG <- readRDS(file.path(OUT_S5, "agg.rds")); M <- AG$META
X <- M[set == "own" & method %in% PRIM8, .(dgp, scen, N, ntime, arm, method, j, d)]
X[, `:=`(n = AG$A$n[1, j], k = AG$A$k[1, j], kt = AG$A$kt[1, j], nd = AG$ND[1, d])]
X[method == "glmer_lrt", kt := NA_real_]                      # the LRT has no t version
P <- fread(file.path(TABDIR, "V4_S6_nestedloop.csv"))
chk <- merge(X[arm == "H0"], P, by = c("dgp", "scen", "N", "ntime", "method"))
stopifnot(nrow(chk) == nrow(X[arm == "H0"]), max(abs(chk$n - chk$n_obs)) == 0, max(abs(chk$nd - chk$n_den)) == 0,
          max(abs(fifelse(chk$n > 0, 100 * chk$k / chk$n, 0) - fifelse(chk$n_obs > 0, chk$rej, 0))) < 1e-9)
fwrite(X[, .(dgp, scen, N, ntime, arm, method, n, nd, k, kt)], file.path(TABDIR, "V4_S6_dp_rates.csv"))
cat(sprintf("%d rows; checked against V4_S6_nestedloop.csv (n, nd, normal rate): identical\n", nrow(X)))
