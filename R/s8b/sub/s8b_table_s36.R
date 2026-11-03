#!/usr/bin/env Rscript
# =============================================================================
# s8b_table_s36.R -- S8b sub-task, plan item DEC-T4 (technical decision): Table S36 of the supplement, the performance
# measures of Morris et al. (2019) that the paper does not tabulate elsewhere, by method x mechanism, main grid (M1-M5)
# and M1c, under the alternative (label mmlb, arm H1), usable set U0, own set:
#   bias / empirical SE, model SE / empirical SE - 1 (%), coverage of the 95% interval (%),
#   each averaged over design points with at least x:min_obs_dp usable replicates.
# MCSE: for each cell, the mean over design points of the per-point MCSE columns of the source (mcse_bias / empse,
#   mcse_relerr, 100 * mcse_cover). Whatever the correlation between design points, sd(mean) <= mean(sd), so this is an
#   upper bound for the MCSE of the cell; the caption gives the maximum of the bound over the cells, per measure.
# Input (read-only): tables/variant_cv/V4_perf_primary.csv (variant with the convergence check of the marginalized
#   model and Kauermann-Carroll from fewclust, D-60), tables/variant_cc/V4_perf_primary.csv, tables/V4_perf_primary.csv
#   (only for the provenance check), tables/V5_S7f_kcfc_gate.csv, tables/V5_S7f_numbers.csv (x:min_obs_dp).
# Outputs: tables/V5_S8b_sub_s36.csv (long table), tables/V5_S8b_sub_s36_keys.csv (keys s36|...),
#          tables/V5_S8b_sub_s36_checks.csv. Run from ~/MMM: Rscript --vanilla R/s8b/sub/s8b_table_s36.R
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
SRC <- "tables/variant_cv/V4_perf_primary.csv"
CHK <- list(); chk <- function(name, ok, detail = "") { CHK[[length(CHK) + 1]] <<- data.table(check = name, pass = isTRUE(ok), detail = detail)
  cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
OUT <- list()
put <- function(key, value, kind, how, mcse = NA_real_, source = SRC) {
  stopifnot(length(value) == 1, !is.na(value))
  OUT[[key]] <<- data.table(key_ = key, value = as.character(value), mcse = if (is.na(mcse)) "" else as.character(mcse), kind = kind, source = source, how = how) }
REG <- fread("tables/V5_S7f_numbers.csv", colClasses = "character")
MIN_OBS <- as.integer(REG[id == "x:min_obs_dp", value]); stopifnot(length(MIN_OBS) == 1, MIN_OBS == 20L)
METH <- c("gee_fg", "gee_md", "mmm", "glmer", "shrink", "gee", "gee_kc")
MECH <- c(GenBinaryY = "GB", `GLMM-native` = "NAT", `GLMM-AR1` = "AR1")
## ---- 1. provenance: this file is the kcfc variant -----------------------------------------------------------------------
P  <- fread(SRC); PC <- fread("tables/variant_cc/V4_perf_primary.csv"); PS <- fread("tables/V4_perf_primary.csv")
K <- c("dgp", "scen", "model", "N", "ntime", "label", "arm", "method", "usab", "set")
VC <- c("n_obs", "rej", "bias", "empse", "modse", "cover", "p_usable")
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-12)
diffcells <- function(A, B) { M <- merge(A[, c(K, VC), with = FALSE], B[, c(K, VC), with = FALSE], by = K, suffixes = c(".a", ".b"))
  d <- Reduce(`|`, lapply(VC, function(v) !same(M[[paste0(v, ".a")]], M[[paste0(v, ".b")]])))
  list(n_common = nrow(M), D = M[d, .N, by = .(method, scen)][order(method, scen)]) }
d1 <- diffcells(PC, P); d2 <- diffcells(PS, PC)
chk("same rows in kcfc and cc variants", nrow(P) == nrow(PC) && d1$n_common == nrow(P), sprintf("%d rows", nrow(P)))
chk("cv differs from cc only in gee_kc rows (fewclust KC, D-60), gee rows (constant outcomes, D-82) and mmm rows of M3-M5 (refit not converged, D-90)", nrow(d1$D) > 0 && all(d1$D$method %in% c("gee_kc", "gee") | (d1$D$method == "mmm" & d1$D$scen %in% c("M3", "M4", "M5"))),   # D-82, D-90
    paste(unique(d1$D$method), collapse = ","))
chk("cc differs from the study row set only in mmm rows of M3-M5 (convergence check, D-32)",
    nrow(d2$D) > 0 && all(d2$D$method == "mmm") && all(d2$D$scen %in% c("M3", "M4", "M5")),
    paste(unique(paste(d2$D$method, d2$D$scen)), collapse = ","))
G <- fread("tables/V5_S7f_kcfc_gate.csv"); chk("kcfc build gate (V5_S7f_kcfc_gate.csv) all pass", all(G$pass), sprintf("%d checks", nrow(G)))
chk("registry tables of the supplement read the cv variant", REG[grepl("^x:tc_T12", id), all(grepl("variant_cv", file))])
## ---- 2. per-point measures ------------------------------------------------------------------------------------------
X <- P[label == "mmlb" & arm == "H1" & usab == "U0" & set == "own" & method %in% METH & scen %in% c(sprintf("M%d", 1:5), "M1c")]
chk("one row per design point and method", X[, .(cnt = .N), by = .(dgp, scen, N, ntime, method)][, all(cnt == 1L)])
X[, grid := fifelse(scen == "M1c", "M1c", "main")]
rr <- X[is.finite(relerr_modse) & is.finite(modse) & is.finite(empse) & empse > 0]
chk("relerr_modse is 100 (modse / empse - 1), in percent", rr[, max(abs(relerr_modse - 100 * (modse / empse - 1)))] < 1e-8,
    sprintf("max abs diff %.2e over %d rows", rr[, max(abs(relerr_modse - 100 * (modse / empse - 1)))], nrow(rr)))
chk("cover is a proportion in [0, 1]", X[is.finite(cover), all(cover >= 0 & cover <= 1)])
chk("GLMM on GenBinaryY: bias and coverage missing in the source", X[method == "glmer" & dgp == "GenBinaryY", all(is.na(bias)) && all(is.na(cover))])
X[, ok := n_obs >= MIN_OBS]
X[, `:=`(rb = bias / empse, mc_rb = mcse_bias / empse, relerr = relerr_modse, mc_relerr = mcse_relerr, cov = 100 * cover, mc_cov = 100 * mcse_cover)]
X[method == "glmer" & dgp == "GenBinaryY", `:=`(rb = NA_real_, mc_rb = NA_real_, cov = NA_real_, mc_cov = NA_real_)]
MS <- c(rb = "mc_rb", relerr = "mc_relerr", cov = "mc_cov")
L <- rbindlist(lapply(names(MS), function(m) X[, {
  u <- ok & is.finite(get(m)) & is.finite(get(MS[[m]]))
  .(measure = m, n_dp = .N, n_used = sum(u), value = mean(get(m)[u]), mcse_bound = mean(get(MS[[m]])[u]), sd_points = sd(get(m)[u]))
  }, by = .(method, dgp, grid)]))
L[n_used == 0L, `:=`(value = NA_real_, mcse_bound = NA_real_)]
L[, mech := MECH[dgp]]
setorder(L, measure, grid, method, mech)
fwrite(L, "tables/V5_S8b_sub_s36.csv")
chk("every cell defined except GLMM on GenBinaryY (bias, coverage)", L[!(method == "glmer" & dgp == "GenBinaryY" & measure != "relerr"), all(n_used > 0)])
chk("design points per grid and mechanism", TRUE, paste(L[measure == "relerr" & method == "mmm", sprintf("%s/%s=%d", grid, mech, n_dp)], collapse = " "))
## ---- 3. keys ----------------------------------------------------------------------------------------------------------
for (i in seq_len(nrow(L))) { r <- L[i]; k <- sprintf("s36|%s|%s|%s|%s", r$method, r$mech, r$grid, r$measure)
  if (is.na(r$value)) put(k, "n.a.", "text", "not defined: the GLMM has no marginal estimand on GenBinaryY")
  else put(k, r$value, "sig2", sprintf("mean over %d of %d design points with n_obs >= %d; mcse = mean of per-point MCSE (upper bound)", r$n_used, r$n_dp, MIN_OBS), mcse = r$mcse_bound) }
for (m in names(MS)) { w <- L[measure == m & is.finite(mcse_bound)][which.max(mcse_bound)]
  put(sprintf("s36|mcse_max|%s", m), w$mcse_bound, "sig2", sprintf("max over the cells of Table S36 of the mean per-point MCSE (at %s %s %s)", w$method, w$mech, w$grid)) }
for (g in c("main", "M1c")) put(sprintf("s36|ndp|%s", g), L[measure == "relerr" & method == "mmm" & grid == g & mech == "GB", n_dp], "int",
                                 "design points per mechanism in the grid (H1)")
for (g in c("main", "M1c")) put(sprintf("s36|nused_min|%s", g), L[n_used > 0 & grid == g, min(n_used)], "int", "smallest number of design points entering a cell of the grid")
put("s36|nused_min", L[n_used > 0, min(n_used)], "int", "smallest number of design points entering a cell")
K36 <- rbindlist(OUT); setnames(K36, "key_", "key"); fwrite(K36, "tables/V5_S8b_sub_s36_keys.csv")
CK <- rbindlist(CHK); fwrite(CK, "tables/V5_S8b_sub_s36_checks.csv")
print(dcast(L, method + grid ~ measure + mech, value.var = "value"), digits = 3)
print(L[, .(max_bound = max(mcse_bound, na.rm = TRUE)), by = measure])
print(L[, .(min_used = min(n_used), max_used = max(n_used), n_dp = max(n_dp)), by = .(grid, mech)])
print(d1$D); print(d2$D)
cat("keys:", nrow(K36), "| checks failed:", sum(!CK$pass), "\n"); stopifnot(all(CK$pass))
