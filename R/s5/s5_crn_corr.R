# =============================================================================
# s5_crn_corr.R -- S5 doubt (b): correlation of the rejection indicators of the SAME replicate
# across cells, by relation between the cells. Explains why the block bootstrap SE of a pooled
# rate equals the SE under independent cells (P3). H0, U0 own, methods mmm and gee_fg, scenarios M1 and M3.
# Output: tables/V4_S5_fu_b_crn_corr.csv
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R"))
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[!is.na(term) & label == "manual" & scen %in% c("M1", "M3") & method %in% c("mmm", "gee_fg"), .(dgp, scen, N, ntime, method, sim_id, p, U0)]
V[, r := fifelse(U0 %in% TRUE, as.numeric(p < 0.05), NA_real_)]
out <- list()
for (sc in c("M1", "M3")) for (mm in c("mmm", "gee_fg")) {
  W <- dcast(V[scen == sc & method == mm], sim_id ~ dgp + N + ntime, value.var = "r")
  key <- tstrsplit(names(W)[-1], "_", fixed = TRUE); K <- data.table(dgp = key[[1]], N = as.integer(key[[2]]), ntime = as.integer(key[[3]]))
  C <- suppressWarnings(cor(as.matrix(W[, -1]), use = "pairwise.complete.obs"))
  ij <- which(upper.tri(C), arr.ind = TRUE)
  R <- data.table(c = C[ij], same_dgp = K$dgp[ij[, 1]] == K$dgp[ij[, 2]], same_N = K$N[ij[, 1]] == K$N[ij[, 2]], same_T = K$ntime[ij[, 1]] == K$ntime[ij[, 2]],
                  dpair = paste(pmin(K$dgp[ij[, 1]], K$dgp[ij[, 2]]), pmax(K$dgp[ij[, 1]], K$dgp[ij[, 2]]), sep = "~"))
  R[, relation := fifelse(same_dgp & same_N & !same_T, "same mechanism, same N, other T",
                  fifelse(same_dgp & !same_N & same_T, "same mechanism, other N, same T",
                  fifelse(same_dgp, "same mechanism, other N and T",
                  fifelse(same_N & same_T, paste0("same N and T, mechanisms ", dpair), "other mechanism, other N or T"))))]
  out[[length(out) + 1]] <- R[, .(scen = sc, method = mm, pairs = .N, mean_corr = mean(c, na.rm = TRUE), q90 = quantile(c, 0.9, na.rm = TRUE)), by = relation]
}
X <- rbindlist(out); print(X[order(scen, method, -mean_corr)]); fwrite(X, file.path(TABDIR, "V4_S5_fu_b_crn_corr.csv"))
