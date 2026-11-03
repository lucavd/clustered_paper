source(path.expand("~/MMM/R/s5/s5_common.R")); source(file.path(REPO, "R/s5/s5_engine.R"))
AG <- readRDS(file.path(OUT_S5, "agg.rds")); META <- AG$META
CB <- readRDS(file.path(OUT_S5, "claims_boot.rds")); PBIN <- readRDS(file.path(OUT_S5, "pbin.rds"))
TO <- function(x) file.path(TABDIR, x)
mk0 <- function() { cl <- META[, .(d, dgp, scen, N, ntime, arm, method, set, tok)]
  cl[, `:=`(n = AG$A$n[1, ], k = AG$A$k[1, ], nd = AG$ND[1, d])]; cl }
cl0 <- mk0()
# ---- (g) LRT usability by N and scenario group, v3 vs v4 --------------------------------------------
U4 <- cl0[scen != "M1c" & arm == "H0" & set == "own" & method %in% PRIM8][, pu := 100 * n / nd]
REV <- fread(file.path(REPO, "rev2/REV_cells.csv"))[label == "manual" & set == "own"]
att <- REV[, .(nd = max(n_attempt)), by = .(dgp, model, N, ntime)]
g3 <- CJ(dgp = DGPS, model = 1:5, N = sort(unique(REV$N)), ntime = sort(unique(REV$ntime)), method = PRIM8)
g3 <- merge(merge(g3, att, by = c("dgp", "model", "N", "ntime")), REV[, .(dgp, model, N, ntime, method, n_usable)], by = c("dgp", "model", "N", "ntime", "method"), all.x = TRUE)
g3[is.na(n_usable), n_usable := 0][, pu := 100 * n_usable / nd]
grp <- function(s) fifelse(s %in% c(1, 2, "M1", "M2"), "M1-M2", "M3-M5")
LU <- merge(U4[N %in% c(10, 25, 35, 45), .(v4 = mean(pu)), by = .(method, N, ntime_grp = fifelse(ntime <= 3, "T<=3", "T>=4"), sgrp = grp(scen))],
            g3[N %in% c(10, 25, 35, 45), .(v3 = mean(pu)), by = .(method, N, ntime_grp = fifelse(ntime <= 3, "T<=3", "T>=4"), sgrp = grp(model))],
            by = c("method", "N", "ntime_grp", "sgrp"))
print(LU[method %in% c("glmer_lrt", "gee_md", "gee_fg") & N == 35]); fwrite(LU, TO("V4_S5_fu_g_usable.csv"))
VN <- merge(U4[, .(v4_pu = mean(pu), v4_va = mean(100 * (n - k) / nd)), by = .(dgp, method, N)],
            g3[, .(v3_pu = mean(pu)), by = .(dgp, method, N)], by = c("dgp", "method", "N"))
fwrite(VN, TO("V4_S5_fu_g_validity_byN.csv"))
