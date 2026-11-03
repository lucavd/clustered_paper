# =============================================================================
# s2_b_degenerate.R -- verification of the operational definition of D-11:
# a replicate is non-estimable iff its design has a single treatment arm or an
# arm never observed beyond time 0. Regenerates nclust and trt only (first two
# design draws, identical in the three mechanisms) for all 480,000 keys.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
grid <- CJ(model = 1:5, N = N_LEVELS, ntime = T_LEVELS)
one <- function(m, N, Tt) {
  out <- logical(1000)
  for (i in 1:1000) {
    set_rep_seed(m, i)
    nclust <- ifelse(runif(N) < 0.7, Tt, sample(seq_len(Tt - 1), N, replace = TRUE))
    trt <- rbinom(N, 1, 0.5)
    out[i] <- sum(trt == 1) == 0 || sum(trt == 0) == 0 || all(nclust[trt == 1] == 1) || all(nclust[trt == 0] == 1)
  }
  data.table(model = m, N = N, ntime = Tt, sim_id = which(out))
}
DG <- rbindlist(parallel::mcMap(one, grid$model, grid$N, grid$ntime, mc.cores = 40))
lost <- unique(fread(file.path(TABDIR, "V4_S2_missing_ids.csv"))[, .(model, N, ntime, sim_id)])
both <- merge(DG[, deg := TRUE], lost[, lst := TRUE], by = c("model", "N", "ntime", "sim_id"), all = TRUE)
cat(sprintf("degenerate designs %d | lost keys %d | both %d | degenerate not lost %d | lost not degenerate %d\n",
            nrow(DG), nrow(lost), both[deg & lst, .N], both[deg & is.na(lst), .N], both[lst & is.na(deg), .N]))
print(both[is.na(lst) | is.na(deg)])
fwrite(DG, file.path(TABDIR, "V4_S2_degenerate_designs.csv"))
