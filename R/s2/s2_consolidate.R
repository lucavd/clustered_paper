# =============================================================================
# s2_consolidate.R -- S2 summaries read by the report.
# (1) usable shares with every (cell x method) completed (absent = 0 usable);
# (2) exploratory, NOT pre-registered (input for D-04): type I error on H0 cells
#     under U0 (warnings admitted) vs U1 (any warning excluded).
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
CELL <- c("dgp", "model", "N", "ntime", "label")
G <- fread(file.path(TABDIR, "V4_S2_replicates_by_cell.csv"))
U <- fread(file.path(TABDIR, "V4_S2_usable_by_cell.csv"))
METH <- c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt")
full <- G[, .(method = METH), by = c(CELL, "n_gen")]
UC <- merge(full, U[, c(CELL, "method", "n_U0", "n_U1", "n_U2", "n_warn"), with = FALSE], by = c(CELL, "method"), all.x = TRUE)
for (v in c("n_U0", "n_U1", "n_U2", "n_warn")) UC[is.na(get(v)), (v) := 0L]
cat(sprintf("completed rows %d (= 2880 x 8: %s); rows added with zero %d\n", nrow(UC), nrow(UC) == 23040, nrow(UC) - nrow(U)))
P <- UC[, .(n_gen = sum(n_gen), U0 = sum(n_U0) / sum(n_gen), U1 = sum(n_U1) / sum(n_gen), U2 = sum(n_U2) / sum(n_gen),
            warn_among_U0 = sum(n_warn) / sum(n_U0),
            cell_med_warn = median(fifelse(n_U0 > 0, n_warn / n_U0, NA_real_), na.rm = TRUE),
            cells_warn_gt10 = sum(n_U0 > 0 & n_warn / n_U0 > 0.10), cells_warn_gt50 = sum(n_U0 > 0 & n_warn / n_U0 > 0.50),
            cells_zero_U0 = sum(n_U0 == 0)), by = method]
fwrite(P, file.path(TABDIR, "V4_S2_usable_pooled.csv")); print(P, digits = 4)
PT <- UC[, .(U0 = sum(n_U0) / sum(n_gen), U1 = sum(n_U1) / sum(n_gen)), by = .(method, ntime)]
fwrite(PT, file.path(TABDIR, "V4_S2_usable_by_T.csv"))

# ---- (2) exploratory: type I error under U0 vs U1 (H0 cells) ----------------------
A <- readRDS(file.path(OUTDIR, "rows_timetrt.rds"))[label == "manual" & method %in% METH & term == "time:trt"]
WALD <- setdiff(METH, "glmer_lrt")
A[, U0 := fifelse(method %in% WALD, converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p),
                  converged & is.finite(p))]
A[is.na(U0), U0 := FALSE]; A[, U1 := U0 & is.na(status)]
T1 <- A[, .(n0 = sum(U0), k0 = sum(U0 & p < 0.05), n1 = sum(U1), k1 = sum(U1 & p < 0.05)), by = .(dgp, model, N, ntime, method)]
T1[, `:=`(t1_U0 = k0 / n0, t1_U1 = fifelse(n1 > 0, k1 / n1, NA_real_))]
T1[, d := t1_U1 - t1_U0]
T1[, mcse_U1 := sqrt(pmax(t1_U1 * (1 - t1_U1), 1e-12) / pmax(n1, 1))]
fwrite(T1, file.path(TABDIR, "V4_S2_typeI_U0_vs_U1_explore.csv"))
S <- T1[n0 >= 100 & n1 >= 100, .(cells = .N, mean_t1_U0 = mean(t1_U0), mean_t1_U1 = mean(t1_U1), mean_d = mean(d),
                                  q05_d = quantile(d, 0.05), q95_d = quantile(d, 0.95), cells_absd_gt_2mcse = sum(abs(d) > 2 * mcse_U1)), by = method]
fwrite(S, file.path(TABDIR, "V4_S2_typeI_U0_vs_U1_summary.csv")); print(S, digits = 3)
cat(sprintf("cells with n1 < 100 but n0 >= 100, by method:\n")); print(T1[n0 >= 100 & n1 < 100, .N, by = method])

# ---- AR1 Q3.5a detail ---------------------------------------------------------------
AS <- fread(file.path(TABDIR, "V4_S2_truth_ar1_conditional.csv"))
cat("Q3.5a cells with sd < 1e-8, by model x ntime:\n"); print(dcast(AS[, .(k = sum(sd < 1e-8, na.rm = TRUE)), by = .(model, ntime)], model ~ ntime, value.var = "k"))
cat(sprintf("Q3.5a max sd %.4g; max sd / |mean| among H1 %.4g\n", max(AS$sd, na.rm = TRUE), AS[label == "mmlb", max(sd / abs(mean), na.rm = TRUE)]))
