
# s10s_ci_agq.R -- S10 "supplemento": inverted sign-flip interval (Table S27, Section S7) and adaptive quadrature
# (Table S34, Section S12), from the per-replicate rows copied by s10s_build.R. Definitions from the rendered text.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2)
D <- "/home/user/MMM/results_v4/S10/supplemento"; K <- c("dgp","scen","N","ntime","label","sim_id"); CELL <- c("dgp","scen","N","ntime","label")
# ---------------- inverted sign-flip interval
C <- readRDS(file.path(D, "ci_s7d2.rds"))
C[, estd := generated %in% TRUE & !(design_degenerate %in% TRUE)]
C[, const := estd & (sum_y == 0 | sum_y == n_obs)]
C[, np := fifelse(scen %in% c("M2","M5"), 5L, 4L)]
C[, fin := is.finite(lo) & is.finite(hi)]
C[, gu := ar1_gee_ok %in% TRUE & is.finite(ar1_est) & is.finite(ar1_se_fg) & abs(ar1_est) < 10 & ar1_se_fg > 0 & ar1_se_fg < 10]
C[, `:=`(gu_nc = gu & !const, fin_nc = fin & !const)]
C[, gw := 2 * qt(0.975, N - np) * ar1_se_fg]; C[, fw := hi - lo]
C[, Ng := fcase(N == 10, "10", N <= 25, "15-25", N <= 65, "35-65", default = "75-145")]
T27 <- C[label == "mmlb" & estd == TRUE, .(datasets = .N, const = sum(const), finite = 100 * mean(fin), finite_nc = 100 * mean(fin_nc),
            gee = 100 * mean(gu), gee_nc = 100 * mean(gu_nc), finite_amongnc = 100 * mean(fin[!const])), by = .(ntime, Ng)]
fwrite(T27[order(ntime, Ng)], file.path(D, "s10s_ci_T27.csv"))
small <- C[label == "mmlb" & estd == TRUE & ntime <= 3 & N <= 25]
both <- C[label == "mmlb" & estd == TRUE & fin & gu_nc]
txt <- data.table(item = c("PG6_finite_small_H1", "PG6_finite_small_H1_nc", "ci_width_gm_pct", "ci_width_gm_pct_N75",
                          "ci_width_gm_pct_inclconst", "n_both"),
                  value = c(100 * mean(small$fin), 100 * mean(small$fin_nc), 100 * (exp(mean(log(both$fw / both$gw))) - 1),
                            both[N >= 75, 100 * (exp(mean(log(fw / gw))) - 1)],
                            C[label == "mmlb" & estd == TRUE & fin & gu, 100 * (exp(mean(log(fw / gw))) - 1)], nrow(both)))
# binomial MCSE cap of Table S27 (finite share), per row
T27[, mcse_fin := 100 * sqrt(finite / 100 * (1 - finite / 100) / datasets)]
txt <- rbind(txt, data.table(item = c("cap_S27_finite", "cap_S27_gee"), value = c(max(T27$mcse_fin), max(100 * sqrt(T27$gee / 100 * (1 - T27$gee / 100) / T27$datasets)))))
fwrite(txt, file.path(D, "s10s_ci_text.csv")); rm(C); gc()
# ---------------- adaptive quadrature
Rg <- readRDS(file.path(D, "recs_s7g.rds")); Wg <- readRDS(file.path(D, "rows_s7g.rds"))
Rg[, estd := generated %in% TRUE & !(design_degenerate %in% TRUE)][, const := estd & (sum_y == 0 | sum_y == n_obs)]
W <- merge(Wg, Rg[, c(K, "estd", "const"), with = FALSE], by = K)
W <- W[estd == TRUE]
W[, np := fifelse(scen %in% c("M2","M5"), 5L, 4L)]
W[, use := converged %in% TRUE & is.finite(estimate) & is.finite(se) & is.finite(p) & abs(estimate) < 10 & se > 0 & se < 10 & !const]
W[, rt := abs(estimate / se) > qt(0.975, N - np)]
A <- W[, .(n_est = .N, n_use = sum(use), k = sum(use & rt)), by = c(CELL, "nAGQ")]
fwrite(A, file.path(D, "s10s_agq_cells.csv"))
# decisions at 25 nodes vs 10 on the first 100 replicates
d <- dcast(W[sim_id <= 100], dgp + scen + N + ntime + label + sim_id ~ nAGQ, value.var = c("use", "rt"))
fwrite(data.table(item = c("rows_nagq1", "n_both_10_25", "diff_dec_10_25"),
                  value = c(W[nAGQ == 1, .N], d[use_10 & use_25, .N], d[use_10 & use_25 & rt_10 != rt_25, .N])), file.path(D, "s10s_agq_text.csv"))
cat("CI_AGQ OK\n")
