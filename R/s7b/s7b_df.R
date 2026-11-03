# =============================================================================
# s7b_df.R -- S7b part 1: Wald tests with the t at N - 2 degrees of freedom (Li and Redden) next to N - p,
# on the stored estimates and standard errors of all 1,896 design points (1,512 of the study with the MMM
# convergence check, results_v4/variant_cc/S4/rows_v4.rds; 384 of M6-M7, results_v4/S6b/ns/rows_ns.rds).
# Gates GC2 (reproduction of the stored N - p counts and of the stored normal p), GC3 (completeness).
# Outputs: tables/V4_S7b_dp_rates.csv   per design point x arm x method: n, nd, k (normal), kt (N - p), kt2 (N - 2), empse
#          tables/V4_S7b_rule_dp.csv    per design point x strategy (method|n, |t, |t2): rate, mc, flags, pu, upow
#          tables/V4_S7b_rule_summary.csv  by strategy x source (studio / nuovi / tutti) and the rule's choice
#          tables/V4_S7b_df_shift.csv   N - 2 minus N - p rejection rate, mean over design points, block bootstrap MCSE
#          tables/V4_S7b_gates.csv      GC2, GC3
#          results_v4/S7b/rows_h0_marg.rds  per-replicate estimates of H0 GLMM-native marginal methods (for s7b_null.R)
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R"))
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
source(path.expand("~/MMM/R/s11/s11_cv.R"))   # D-90 (S11, D-91): no result for the MMM with the check when its refit does not converge
t0 <- Sys.time()
KEEP <- c("dgp", "scen", "N", "ntime", "label", "method", "sim_id", "term", "estimate", "se", "p", "U0")
A <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))[method %in% PRIM8, ..KEEP]
A[, src := "studio"]
B <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds"))
B <- B[method %in% PRIM8, c(KEEP, "design_degenerate"), with = FALSE][, src := "nuovi"]
cat(sprintf("rows: study %d, new scenarios %d | load %.1f min\n", nrow(A), nrow(B), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
# degenerate designs: study from the S2 table (model, N, ntime, sim_id; M1c has the model and seeds of M1), M6-M7 from the row flag
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
A[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
A[, deg := FALSE]; A[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
B[, deg := design_degenerate %in% TRUE]; B[, design_degenerate := NULL]; B[, model := as.integer(sub("^M", "", scen))]
V <- rbindlist(list(A, B), use.names = TRUE); rm(A, B); invisible(gc(verbose = FALSE))
V <- V[!is.na(term)]
V[, u0 := U0 & !deg]; V[, arm := fifelse(label == "mmlb", "H1", "H0")]
V[, u0_pre := u0]; s9_drop_const(V, "u0", what = "study and M6-M7 rows"); D82_TOUCHED <- V[u0_pre & !u0, unique(method)]; V[, u0_pre := NULL]   # D-82
cat("D-82: methods with usable constant-outcome rows:", paste(D82_TOUCHED, collapse = " "), "\n")
V[, u0_pre := u0]; s11_drop_cv(V, "u0", rows = V$method == "mmm", what = "study and M6-M7 rows"); D90_TOUCHED <- V[u0_pre & !u0, unique(method)]; V[, u0_pre := NULL]   # D-90
V[, np := NP[scen]]
V[, z := estimate / se]
V[, r  := u0 & fifelse(is.na(p), FALSE, p < 0.05)]
V[, rt := u0 & method != LRT & fifelse(is.finite(z), 2 * pt(-abs(z), df = N - np) < 0.05, FALSE)]
V[, rt2 := u0 & method != LRT & fifelse(is.finite(z), 2 * pt(-abs(z), df = N - 2L) < 0.05, FALSE)]
# normal p recomputed from estimate and SE (part of GC2)
pn <- V[u0 & method %in% WALD7 & is.finite(z), .(n = .N, n_bad = sum(abs(2 * pnorm(-abs(z)) - p) > 1e-9), max_dev = max(abs(2 * pnorm(-abs(z)) - p)))]
print(pn)
# ---- per design point ----------------------------------------------------------------------------
ND_A <- unique(V[src == "studio", .(dgp, scen, model, N, ntime)])
ND_A <- merge(ND_A, DEG[, .(ndeg = .N), by = .(model, N, ntime)], by = c("model", "N", "ntime"), all.x = TRUE)
ND_A[is.na(ndeg), ndeg := 0L][, nd := 1000L - ndeg][, src := "studio"]
ND_B <- unique(fread(file.path(TABDIR, "V4_S6b_ns_dp.csv"))[, .(dgp, scen, N, ntime, nd)])[, src := "nuovi"]
ND <- rbind(ND_A[, .(dgp, scen, N, ntime, nd, src)], ND_B)
stopifnot(nrow(ND) == 1896L)
X <- V[, .(n = sum(u0), k = sum(r), kt = sum(rt), kt2 = sum(rt2), empse = sd(estimate[u0]), mean_est = mean(estimate[u0])), by = .(dgp, scen, N, ntime, arm, method)]
G <- CJ(d = seq_len(nrow(ND)), arm = c("H0", "H1"), method = PRIM8)
G <- cbind(ND[G$d], G[, .(arm, method)])
X <- merge(G, X, by = c("dgp", "scen", "N", "ntime", "arm", "method"), all.x = TRUE)
for (v in c("n", "k", "kt", "kt2")) set(X, which(is.na(X[[v]])), v, 0)
X[method == LRT, `:=`(kt = NA_real_, kt2 = NA_real_)]
setorder(X, src, dgp, scen, N, ntime, arm, method)
fwrite(X, file.path(TABDIR, "V4_S7b_dp_rates.csv"))
# ---- strategies and flags ---------------------------------------------------------------------------
L <- rbind(X[, .(dgp, scen, N, ntime, src, arm, method, n, nd, kk = k, ref = "n")],
           X[method != LRT, .(dgp, scen, N, ntime, src, arm, method, n, nd, kk = kt, ref = "t")],
           X[method != LRT, .(dgp, scen, N, ntime, src, arm, method, n, nd, kk = kt2, ref = "t2")])
L[, strategy := paste(method, ref, sep = "|")]
H0 <- flags(L[arm == "H0"])
H1 <- L[arm == "H1", .(dgp, scen, N, ntime, strategy, upow = 100 * kk / nd)]
D <- merge(H0, H1, by = c("dgp", "scen", "N", "ntime", "strategy"))
D[, band := fifelse(ntime <= 3, "T2-3", fifelse(ntime <= 6, "T4-6", "T8-10"))]
D <- D[, .(dgp, scen, N, ntime, band, src, strategy, method, ref, n, nd, rate, mc, assess, hi, lo, pu, upow)]
fwrite(D, file.path(TABDIR, "V4_S7b_rule_dp.csv"))
# ---- GC2: reproduction of the stored tables (n, nd, rate, hi, lo) for ref n and t --------------------------
S1 <- fread(file.path(REPO, "tables/variant_cc/V4_S6_rule_safety_dp.csv"))[, src := "studio"]
S2 <- fread(file.path(TABDIR, "V4_S6b_ns_rule_safety_dp.csv"))[!startsWith(strategy, "mmm_nocc")][, src := "nuovi"]
S <- rbind(S1, S2)
M <- merge(D[ref %in% c("n", "t")], S, by = c("dgp", "scen", "N", "ntime", "strategy", "src"), suffixes = c("", ".s"))
gc2_rows <- nrow(M) == nrow(S) && nrow(M) == nrow(D[ref %in% c("n", "t")])
gc2_eq <- M[!method %in% c(D82_TOUCHED, D90_TOUCHED), all(n == n.s & nd == nd.s & hi == hi.s & lo == lo.s & assess == assess.s &
                  ((is.na(rate) & is.na(rate.s)) | abs(rate - rate.s) < 1e-9) & abs(pu - pu.s) < 1e-9 & abs(upow - upow.s) < 1e-9)]
gc2 <- gc2_rows && gc2_eq && pn$n_bad == 0L   # D-82: reproduction of the stored tables checked on the methods the rule does not change
# ---- GC3: completeness ----------------------------------------------------------------------------------
cmp <- D[ref == "t2", .N, by = strategy]
lost <- V[u0 & method %in% WALD7, sum(!is.finite(z))]
gc3 <- nrow(cmp) == 7L && all(cmp$N == 1896L) && lost == 0L
GT <- data.table(gate = c("GC2", "GC3"), pass = c(gc2, gc3),
  detail = c(sprintf("%d design point x strategy rows compared with the stored tables (%d study, %d new scenarios): n, nd, rate, hi, lo, pu, upow identical: %s; normal p recomputed from estimate and SE on %d usable Wald rows, %d beyond 1e-9 (max deviation %.1e)",
                     nrow(M), nrow(S1), nrow(S2), gc2_rows && gc2_eq, pn$n, pn$n_bad, pn$max_dev),
             sprintf("%d Wald strategies x %s design points with the t at N - 2; usable rows with a non-finite statistic: %d", nrow(cmp), paste(unique(cmp$N), collapse = "/"), lost)))
fwrite(GT, file.path(TABDIR, "V4_S7b_gates.csv")); print(GT)
# ---- summary by strategy and the rule --------------------------------------------------------------------
sm <- function(Z, by) Z[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_na = sum(!assess), max_rate = max(rate[assess], na.rm = TRUE),
                           pu = mean(pu), upow = mean(upow)), by = by]
SM <- rbind(sm(D, c("strategy", "method", "ref", "src")), sm(D, c("strategy", "method", "ref"))[, src := "tutti"],
            sm(D[N == 10], c("strategy", "method", "ref"))[, src := "tutti N=10"],
            sm(D[N <= 25], c("strategy", "method", "ref"))[, src := "tutti N<=25"], fill = TRUE)
SM[, eligible := n_hi == 0]
SM[src == "tutti", recommended := eligible & pu == max(pu[eligible]), by = ref]
setorder(SM, ref, src, -eligible, -pu); fwrite(SM, file.path(TABDIR, "V4_S7b_rule_summary.csv"))
print(SM[src == "tutti"], digits = 3)
# ---- hi flags that differ between t and t2, listed ------------------------------------------------------------
W <- dcast(D[ref %in% c("t", "t2")], dgp + scen + N + ntime + src + method ~ ref, value.var = c("hi", "lo", "rate"))
FL <- W[hi_t != hi_t2 | lo_t != lo_t2]
fwrite(FL, file.path(TABDIR, "V4_S7b_df_flips.csv")); cat("design points whose hi or lo flag changes between N - p and N - 2:", nrow(FL), "\n")
# ---- shift N - 2 minus N - p: mean over design points, block bootstrap by (model, sim_id) --------------------
BI <- c(readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds")), readRDS(file.path(REPO, "results_v4/S6b/ns/boot_idx_67.rds")))
names(BI) <- as.character(1:7)
CNT <- lapply(BI, function(Mx) { Cm <- matrix(0, 2000L, 1000L); for (b in 1:2000) Cm[b, ] <- tabulate(Mx[b, ], 1000L); rbind(1, Cm) })
H <- V[arm == "H0" & method %in% WALD7 & u0 == TRUE]
H[, cell := paste(dgp, scen, N, ntime, method)]
shift_boot <- function(Hs) {
  out <- list()
  for (m in sort(unique(Hs$model))) {
    h <- Hs[model == m]; cl <- unique(h$cell)
    U <- matrix(0, length(cl), 1000L); Dm <- U; ix <- cbind(match(h$cell, cl), h$sim_id)
    U[ix] <- 1; Dm[ix] <- as.numeric(h$rt2) - as.numeric(h$rt)
    Wm <- t(CNT[[as.character(m)]])              # 1000 x 2001
    out[[length(out) + 1L]] <- data.table(cell = cl, method = sub("^.* ", "", cl), N = as.integer(sub("^\\S+ \\S+ (\\d+) .*$", "\\1", cl)),
                                          rate = 100 * (Dm %*% Wm) / (U %*% Wm))
  }
  rbindlist(out)
}
R <- shift_boot(H)
rcols <- grep("^rate", names(R), value = TRUE)
agg <- function(Rs, lab) Rs[, { M <- as.matrix(.SD); mm <- colMeans(M, na.rm = TRUE); .(group = lab, n_dp = .N, shift = mm[1], mcse = sd(mm[-1])) }, by = method, .SDcols = rcols]
SH <- rbind(agg(R, "all design points"), agg(R[N == 10], "N = 10"), agg(R[N <= 25], "N <= 25"), agg(R[N >= 35], "N >= 35"))
fwrite(SH, file.path(TABDIR, "V4_S7b_df_shift.csv")); print(SH, digits = 3)
# ---- rows for part 2 ------------------------------------------------------------------------------------------
saveRDS(V[arm == "H0" & dgp == "GLMM-native" & method %in% MARGINAL_METHODS & u0 == TRUE, .(dgp, scen, model, N, ntime, method, sim_id, estimate, se, src)],
        file.path(OUT_S7B, "rows_h0_marg.rds"), compress = FALSE)
cat(sprintf("done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
