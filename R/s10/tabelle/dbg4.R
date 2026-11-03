
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R")
ALL <- readRDS(file.path(OUT, "s4_all_dp.rds"))
d <- ALL[, .(cnt = .N), by = c(KD, "strategy")][cnt > 1L]; print(d[, .N, by = .(strategy, dgp, scen)])
S2 <- readRDS(file.path(OUT, "s2_study_dp.rds")); print(S2[, .N, by = strategy])
cons <- ALL[, .(nv = uniqueN(paste(nd, nconst))), by = KD][nv > 1L]; print(cons[, .N, by = .(dgp, scen)])
x <- merge(ALL[cons, on = KD][strategy %in% c("flip", "mmm_cc", "gee_exch_fg")], cons[, ..KD], by = KD)
print(head(dcast(x, dgp + scen + N + ntime + label ~ strategy, value.var = c("nd", "nconst")), 12))
