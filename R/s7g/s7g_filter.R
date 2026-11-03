# =============================================================================
# s7g_filter.R -- S7g, doubt 3, prediction PG4: sensitivity of Table 4 to the usability limits (|estimate| < L, 0 < SE < L)
# with L = 5, 10 (as published) and 20, for the 17 strategies of Table 4 that give an estimate and an interval.
# Per-replicate sources and use rules as the code that built Table 4:
#   GEE strategies: results_v4/S7d/rows (R/s7d/s7d_agg.R) and results_v4/S7e/rows $gee (R/s7e/s7e_agg.R): !degenerate, gee_ok,
#     finite estimate and SE, limits.
#   Marginalized model with the check, GLMM Wald, shrinkage GLM: results_v4/variant_cc/S4/rows_v4.rds, results_v4/S6b/ns/rows_ns.rds
#     (u0_flag: converged, finite estimate, SE and p, limits; not degenerate); Copula-AR1: $study and $cc of the S7e rows.
# H0 arm only; t reference with N - p degrees of freedom; exact Clopper-Pearson flags (D-63), assessed with >= 20 usable;
# usable share = usable / replicates with an estimable design. Check: with L = 10 the clearly-above counts and mean usable
# shares equal tables/V5_S7f_table4.csv for every strategy.
# Outputs: tables/V5_S7g_filter_dp.csv, tables/V5_S7g_filter_t4.csv, tables/V5_S7g_filter_check.csv.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(8L); setwd("/home/user/MMM"); t0 <- Sys.time()
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
source(path.expand("~/MMM/R/s11/s11_cv.R"))   # D-90 (S11, D-91): no result for the MMM with the check when its refit does not converge
cp_lo <- function(k, n) fifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
LIM <- c(5, 10, 20); K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
lim_ok <- function(e, s, L) is.finite(e) & is.finite(s) & abs(e) < L & s > 0 & s < L
add <- function(X, s) rbindlist(lapply(LIM, function(L) {
  u <- X$ok & lim_ok(X$estimate, X$se, L); r <- u & 2 * pt(-abs(X$estimate / X$se), X$N - X$np) < 0.05
  data.table(X[, .(dgp, scen, N, ntime)], u = u, r = r)[, .(n = sum(u), k = sum(r, na.rm = TRUE)), by = .(dgp, scen, N, ntime)][, `:=`(strategy = s, L = L)] }))
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
OUT <- list()
## ---- S7d rows (H0 files) ----
f7d <- list.files("results_v4/S7d/rows", pattern = "_manual_c[0-9]+\\.rds$", full.names = TRUE)
L7 <- parallel::mclapply(f7d, readRDS, mc.cores = 12L)
RECS <- rbindlist(lapply(L7, `[[`, "recs"), fill = TRUE); ROWS <- rbindlist(lapply(L7, `[[`, "rows"), fill = TRUE); rm(L7); invisible(gc(verbose = FALSE))
cat(sprintf("S7d H0: %d files, %d recs, %d rows (%.1f min)\n", length(f7d), nrow(RECS), nrow(ROWS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
ROWS <- merge(ROWS, RECS[, c(K6, "design_degenerate"), with = FALSE], by = K6); ROWS[, np := NP[scen]]
ROWS[, okc := TRUE]; s9_drop_const(ROWS, "okc", what = "S7d H0 GEE rows")   # D-82
for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) { if (ty == "cr2" && cn == "ar1") next
  X <- ROWS[cor == cn, .(dgp, scen, N, ntime, estimate = est, se = get(SEC[[ty]]), ok = !design_degenerate & gee_ok & okc, np)]
  OUT[[length(OUT) + 1]] <- add(X, paste0("gee_", cn, "_", ty)) }
ND <- RECS[, .(nd = sum(!design_degenerate)), by = .(dgp, scen, N, ntime)]
rm(ROWS); invisible(gc(verbose = FALSE))
## ---- study rows ----
KEEP <- c(K6, "method", "term", "estimate", "se", "p", "converged")
SA <- readRDS("results_v4/variant_cc/S4/rows_v4.rds")[label == "manual" & method %in% c("mmm", "glmer", "shrink") & term %in% "time:trt", ..KEEP]
SB <- readRDS("results_v4/S6b/ns/rows_ns.rds")[label == "manual" & method %in% c("mmm", "glmer", "shrink") & term %in% "time:trt", ..KEEP]
S <- rbind(SA, SB); rm(SA, SB); invisible(gc(verbose = FALSE))
S <- merge(S, RECS[, c(K6, "design_degenerate"), with = FALSE], by = K6); S[, np := NP[scen]]
S[, ok := !design_degenerate & converged %in% TRUE & is.finite(p)]
s9_drop_const(S, "ok", what = "study rows")   # D-82
s11_drop_cv(S, "ok", rows = S$method == "mmm", what = "study rows")   # D-90
for (m in c("mmm", "glmer", "shrink")) OUT[[length(OUT) + 1]] <- add(S[method == m], c(mmm = "mmm_cc", glmer = "glmer", shrink = "shrink")[[m]])
rm(S); invisible(gc(verbose = FALSE))
## ---- Copula-AR1 ----
fe <- list.files("results_v4/S7e/rows", pattern = "_manual_c[0-9]+\\.rds$", full.names = TRUE)
LE <- parallel::mclapply(fe, readRDS, mc.cores = 12L)
RE <- rbindlist(lapply(LE, `[[`, "recs"), fill = TRUE); ST <- rbindlist(lapply(LE, `[[`, "study"), fill = TRUE)
CC <- rbindlist(lapply(LE, `[[`, "cc"), fill = TRUE); GE <- rbindlist(lapply(LE, `[[`, "gee"), fill = TRUE); rm(LE); invisible(gc(verbose = FALSE))
for (D in list(RE, ST, CC, GE)) if (nrow(D)) D[, `:=`(dgp = "Copula-AR1", label = "manual")]
ND <- rbind(ND, RE[, .(nd = sum(!design_degenerate)), by = .(dgp, scen, N, ntime)])
ST <- ST[method %in% c("mmm", "glmer", "shrink") & term %in% "time:trt"]
Sm <- merge(ST[method == "mmm"], CC[A_ok == TRUE & A_has_term == TRUE, c(K6, "cc_est", "cc_se", "cc_p"), with = FALSE], by = K6, all.x = TRUE)
Sm[!is.na(cc_est), `:=`(estimate = cc_est, se = cc_se, p = cc_p)]; Sm[, method := "mmm_cc"]
ST <- rbind(ST[method != "mmm"], Sm[, names(ST), with = FALSE])
ST <- merge(ST, RE[, c(K6, "design_degenerate"), with = FALSE], by = K6); ST[, np := NP[scen]]
ST[, ok := !design_degenerate & converged %in% TRUE & is.finite(p)]
s9_drop_const(ST, "ok", what = "copula study rows")   # D-82
s11_drop_cv(ST, "ok", rows = ST$method == "mmm_cc", what = "copula study rows")   # D-90
for (m in c("mmm_cc", "glmer", "shrink")) OUT[[length(OUT) + 1]] <- add(ST[method == m], m)
GE <- merge(GE, RE[, c(K6, "design_degenerate"), with = FALSE], by = K6); GE[, np := NP[scen]]
GE[, okc := TRUE]; s9_drop_const(GE, "okc", what = "copula GEE rows")   # D-82
for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) { if (ty == "cr2" && cn == "ar1") next
  X <- GE[cor == cn, .(dgp, scen, N, ntime, estimate = est, se = get(SEC[[ty]]), ok = !design_degenerate & gee_ok %in% TRUE & okc, np)]
  OUT[[length(OUT) + 1]] <- add(X, paste0("gee_", cn, "_", ty)) }
## ---- per point and Table 4 per L ----
DP <- rbindlist(OUT)
G <- CJ(i = seq_len(nrow(ND)), strategy = unique(DP$strategy), L = LIM); G <- cbind(ND[G$i], G[, .(strategy, L)])
DP <- merge(G, DP, by = c("dgp", "scen", "N", "ntime", "strategy", "L"), all.x = TRUE); DP[is.na(n), `:=`(n = 0L, k = 0L)]
DP[, `:=`(assess = n >= 20, rate = fifelse(n > 0, 100 * k / n, NA_real_), pu = 100 * n / nd)]
DP[, hi := assess & 100 * cp_lo(k, n) > 7.5]
stopifnot(DP[, .N, by = .(strategy, L)][, all(N == 2376L)])
fwrite(DP, "tables/V5_S7g_filter_dp.csv")
T4 <- DP[, .(hi_all = sum(hi), pu = mean(pu), mean_rate = mean(rate, na.rm = TRUE)), by = .(strategy, L)]
T4[, eligible := hi_all == 0L]; T4[, chosen := eligible & pu == max(pu[eligible]), by = L]
setorder(T4, L, hi_all, -pu); fwrite(T4, "tables/V5_S7g_filter_t4.csv")
P4 <- fread("tables/V5_S7f_table4.csv")
CK <- merge(T4[L == 10], P4[, .(strategy, hi_pub = hi_all, pu_pub = pu)], by = "strategy")
CK[, ok := hi_all == hi_pub & abs(pu - pu_pub) < 1e-9]; fwrite(CK, "tables/V5_S7g_filter_check.csv")
print(CK[, .(strategy, hi_all, hi_pub, pu, pu_pub, ok)], digits = 6)
print(dcast(T4, strategy ~ L, value.var = "hi_all")); print(T4[eligible == TRUE])
cat(sprintf("check with L = 10: %d of %d strategies reproduce Table 4; %.1f min\n", CK[, sum(ok)], nrow(CK), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
