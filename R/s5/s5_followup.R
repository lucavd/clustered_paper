# =============================================================================
# s5_followup.R -- S5: investigation of the doubts raised by the first pass (report sec. 3).
# (a) G2 mismatch on mcse_power_sa; (b) why the block SE equals the independent SE (P3);
# (c) bootstrap bias of non-smooth statistics; (d) tie rule and cells won (C04);
# (f) RMSE averaged on own vs common cells (C19a); (g) LRT usability by N, v3 vs v4 (C23).
# Outputs: tables/V4_S5_fu_*.csv
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R")); source(file.path(REPO, "R/s5/s5_engine.R"))
AG <- readRDS(file.path(OUT_S5, "agg.rds")); META <- AG$META
CB <- readRDS(file.path(OUT_S5, "claims_boot.rds")); PBIN <- readRDS(file.path(OUT_S5, "pbin.rds"))
TO <- function(x) file.path(TABDIR, x)
# ---- (a) ------------------------------------------------------------------------------------
P <- fread(file.path(TABDIR, "V4_perf_all.csv.gz"))[usab == "U0" & set == "own" & label == "mmlb"]
M <- META[arm == "H1" & set == "own"]
kc <- match(M$j, as.integer(colnames(AG$KSA)))
M[, se_sa := apply(AG$KSA[-1, kc, drop = FALSE] / AG$A$n[-1, j, drop = FALSE], 2, sd, na.rm = TRUE)]
M[, n_fin := apply(AG$KSA[-1, kc, drop = FALSE] / AG$A$n[-1, j, drop = FALSE], 2, function(x) sum(is.finite(x)))]
X <- merge(M, P[, .(dgp, scen, N, ntime, method, mcse_power_sa, power_sa, n_obs)], by = c("dgp", "scen", "N", "ntime", "method"))
bad <- X[xor(is.finite(se_sa), is.finite(mcse_power_sa)) | (is.finite(se_sa) & abs(se_sa - mcse_power_sa) > 1e-12)]
print(bad[, .(dgp, scen, N, ntime, method, n_obs, power_sa, se_sa, mcse_power_sa, n_fin)])
fwrite(bad[, .(dgp, scen, N, ntime, method, n_obs, power_sa, se_sa, mcse_power_sa, n_fin)], TO("V4_S5_fu_a_g2.csv"))
# ---- (c) bootstrap bias and basic interval ------------------------------------------------------
bias <- rbindlist(lapply(c("main", "boby"), function(v) { MAT <- CB[[v]]
  nm <- grep("^C[0-9]", colnames(MAT), value = TRUE)
  rbindlist(lapply(nm, function(z) { x <- unname(MAT[-1, z]); e <- unname(MAT[1, z])
    q <- quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    data.table(variant = v, name = z, est = e, boot_mean = mean(x, na.rm = TRUE), bias = mean(x, na.rm = TRUE) - e,
               se = sd(x, na.rm = TRUE), pct_lo = q[1], pct_hi = q[2], basic_lo = 2 * e - q[2], basic_hi = 2 * e - q[1],
               est_inside_pct = e >= q[1] & e <= q[2]) })) }))
fwrite(bias, TO("V4_S5_fu_c_bias.csv"))
print(bias[variant == "main" & est_inside_pct == FALSE, .(name, est, boot_mean, pct_lo, pct_hi, basic_lo, basic_hi)])
# ---- (d) tie rule ----------------------------------------------------------------------------------
mk0 <- function() { cl <- META[, .(d, dgp, scen, N, ntime, arm, method, set, tok)]
  cl[, `:=`(n = AG$A$n[1, ], k = AG$A$k[1, ], kt = AG$A$kt[1, ], kc = AG$A$kc[1, ], e2 = AG$A$e2[1, ], ksa = NA_real_, nd = AG$ND[1, d])]; cl }
cl0 <- mk0()
ALPHA <- sort(PRIM8)                                  # order of the methods in O_level_all.csv
cl0a <- copy(cl0)[, mo := match(method, ALPHA)]; setorder(cl0a, d, arm, set, mo); cl0a[, mo := NULL]
w_v4_first <- engine(cl0a, tie = "first", PBIN = PBIN); w_v4_share <- CB$main[1, ]
OL <- fread(file.path(TABDIR, "O_level_all.csv"))[set == "own"]
OL[, dev := abs(100 * k_rej / n - 5)]
wv3 <- function(D, tie) { if (tie == "share") D[, { mn <- min(dev); w <- method[dev == mn]; .(method = w, w = 1 / length(w)) }, by = .(dgp, model, N, ntime)][, .(w = sum(w)), by = method]
                          else D[, .(method = method[which.min(dev)], w = 1), by = .(dgp, model, N, ntime)][, .(w = sum(w)), by = method] }
v3f <- wv3(OL, "first"); v3s <- wv3(OL, "share")
ties <- OL[, .(nt = sum(dev == min(dev)), who = paste(sort(method[dev == min(dev)]), collapse = "+")), by = .(dgp, model, N, ntime)]
D4 <- cl0[scen != "M1c" & arm == "H0" & set == "own" & method %in% PRIM8 & n > 0][, dev := abs(100 * k / n - 5)]
ties4 <- D4[, .(nt = sum(dev == min(dev)), who = paste(sort(method[dev == min(dev)]), collapse = "+")), by = .(dgp, scen, N, ntime)]
WT <- rbindlist(lapply(PRIM8, function(mm) data.table(method = mm, v3_first = v3f[method == mm, sum(w)], v3_share = v3s[method == mm, sum(w)],
                                                     v4_first = w_v4_first[[paste0("wins|", mm)]], v4_share = w_v4_share[[paste0("wins|", mm)]])))
print(WT); fwrite(WT, TO("V4_S5_fu_d_wins.csv"))
cat("v3 design points with ties:", ties[nt > 1, .N], "| v4:", ties4[nt > 1, .N], "\n")
print(head(ties4[nt > 1, .N, by = who][order(-N)], 8)); print(head(ties[nt > 1, .N, by = who][order(-N)], 8))
fwrite(rbind(ties[nt > 1, .N, by = who][, version := "v3"], ties4[nt > 1, .N, by = who][, version := "v4"]), TO("V4_S5_fu_d_ties.csv"))
# ---- (f) RMSE on own cells vs common cells ----------------------------------------------------------
H1 <- cl0[scen != "M1c" & arm == "H1" & set == "own" & method %in% WALD7 & tok == TRUE & n >= 2][, rmse := sqrt(e2 / n)]
own <- H1[, .(rmse_own = mean(rmse), cells_own = .N), by = .(dgp, method)]
cm <- H1[, .(m = uniqueN(method)), by = .(dgp, scen, N, ntime)]
H1 <- merge(H1, H1[, .(mx = uniqueN(method)), by = dgp], by = "dgp")
com <- merge(H1, cm, by = c("dgp", "scen", "N", "ntime"))[m == mx, .(rmse_common = mean(rmse), cells_common = .N), by = .(dgp, method)]
RM <- merge(own, com, by = c("dgp", "method")); print(RM[order(dgp, rmse_common)]); fwrite(RM, TO("V4_S5_fu_f_rmse.csv"))
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
