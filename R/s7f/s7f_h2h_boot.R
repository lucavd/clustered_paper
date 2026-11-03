#!/usr/bin/env Rscript
# s7f_h2h_boot.R -- S7f (second adversarial review): the direct comparison of the two strategies of the rule with the
# block bootstrap declared in Section 2.5 and the relevance margins of Section 2.5. Per-replicate reconstruction adapted
# from ~/MMM_rev2_v5/rev2_v5_recompute.R (Revisore n.2, 2026-10-26). Inputs read-only (results_v4 rows; exact-interval
# per-point files for the gate). Outputs: tables/V5_S7f_h2h_boot.csv, tables/V5_S7f_T2_ar1_exch_identity.csv,
# results_v4/S7f/per_replicate_gee_ar1_fg_mmm_cc.rds. Run from ~/MMM after s7f_cpflags.R.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L)
REPO <- "/home/user/MMM"; OUT <- "/home/user/MMM/results_v4/S7f"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE); setwd(REPO)
source(path.expand("~/MMM/R/s9/s9_const.R"))   # D-82 (S9, D-83): constant-outcome datasets are not usable
source(path.expand("~/MMM/R/s11/s11_cv.R"))   # D-90 (S11, D-91): no result for the MMM with the check when its refit does not converge
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L, M3 = 4L, M4 = 4L, M5 = 5L, M6 = 4L, M7 = 4L)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
rd <- function(dir, gee_name) {
  f <- list.files(file.path(REPO, dir), pattern = "\\.rds$", full.names = TRUE)
  L <- parallel::mclapply(f, function(x) { y <- readRDS(x); g <- y[[gee_name]][cor %in% c("ar1", "exch")]
       list(recs = y$recs[, c(K6, "design_degenerate", "p_flip"), with = FALSE],
            gee = g[, c(K6, "cor", "gee_ok", "est", "se_fg"), with = FALSE],
            study = if (!is.null(y$study)) y$study else NULL, cc = if (!is.null(y$cc)) y$cc else NULL) }, mc.cores = 8L)
  list(recs = rbindlist(lapply(L, `[[`, "recs")), gee = rbindlist(lapply(L, `[[`, "gee")),
       study = rbindlist(lapply(L, `[[`, "study"), fill = TRUE), cc = rbindlist(lapply(L, `[[`, "cc"), fill = TRUE))
}
S7d <- rd("results_v4/S7d/rows", "rows"); tm("S7d read")
S7e <- rd("results_v4/S7e/rows", "gee"); tm("S7e read")
RECS <- rbind(S7d$recs, S7e$recs); GEE <- rbind(S7d$gee, S7e$gee)
cat("recs", nrow(RECS), "unique", uniqueN(RECS[, ..K6]), " gee rows", nrow(GEE), "\n")
GEE <- merge(GEE, RECS[, c(K6, "design_degenerate"), with = FALSE], by = K6)
GEE[, np := NP[scen]]
# usability as written in Section 2.5 of the manuscript
GEE[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(est) & abs(est) < 10 & is.finite(se_fg) & se_fg > 0 & se_fg < 10]
s9_drop_const(GEE, "use", what = "GEE AR(1) and exchangeable FG")   # D-82
GEE[, z := est / se_fg]; GEE[, rej := use & 2 * pt(-abs(z), N - np) < 0.05]
# T = 2: AR(1) vs exchangeable identity (one-parameter correlations coincide with two occasions)
W2 <- dcast(GEE[ntime == 2], dgp + scen + N + label + sim_id ~ cor, value.var = c("est", "se_fg", "use", "rej"))
id2 <- W2[, .(n_both = sum(use_ar1 & use_exch, na.rm = TRUE), max_abs_est_diff = max(abs(est_ar1 - est_exch)[use_ar1 & use_exch], na.rm = TRUE),
              max_rel_se_diff = max(abs(se_fg_ar1 / se_fg_exch - 1)[use_ar1 & use_exch], na.rm = TRUE),
              rej_discordant = sum(rej_ar1 != rej_exch, na.rm = TRUE), use_discordant = sum(use_ar1 != use_exch, na.rm = TRUE), n = .N)]
fwrite(id2, "tables/V5_S7f_T2_ar1_exch_identity.csv"); print(id2)
G <- GEE[cor == "ar1", c(K6, "design_degenerate", "use", "rej", "z"), with = FALSE][, strategy := "gee_ar1_fg"]
rm(GEE); invisible(gc())
# ---- MMM with convergence check: study rows (main grid + M1c), S6b (M6-M7), S7e (copula) --------------------------
SA <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))[method == "mmm" & term %in% "time:trt"]
tm(sprintf("variant_cc rows: %d; columns: %s", nrow(SA), paste(names(SA), collapse = ",")))
SB <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds"))[method == "mmm" & term %in% "time:trt"]
tm(sprintf("S6b rows: %d", nrow(SB)))
keep <- c(K6, "estimate", "se", "p", "converged", "U0")
for (nm in c("p", "converged")) { if (!nm %in% names(SA)) SA[, (nm) := NA]; if (!nm %in% names(SB)) SB[, (nm) := NA] }
cat("S6b columns:", paste(names(SB), collapse = ","), "\n")
SAB <- rbind(SA[, ..keep], SB[, ..keep]); rm(SA, SB); invisible(gc())
# copula: mmm rows replaced by the cc refit where present (logic of R/s7e/s7e_agg.R lines 54-55)
ST <- S7e$study[method == "mmm" & term %in% "time:trt"]; CC <- S7e$cc
if (nrow(CC) > 0 && all(c("A_ok", "A_has_term") %in% names(CC))) {
  Sm <- merge(ST, CC[A_ok == TRUE & A_has_term == TRUE, c(K6, "cc_est", "cc_se", "cc_p"), with = FALSE], by = K6, all.x = TRUE)
  Sm[!is.na(cc_est), `:=`(estimate = cc_est, se = cc_se, p = cc_p)]
} else Sm <- ST
cat("copula cc replaced:", if ("cc_est" %in% names(Sm)) Sm[!is.na(cc_est), .N] else 0, "of", nrow(Sm), "\n")
Sm[, U0 := converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p)]; Sm[is.na(U0), U0 := FALSE]
M <- rbind(SAB, Sm[, ..keep])
M <- merge(RECS[, c(K6, "design_degenerate"), with = FALSE], M, by = K6, all.x = TRUE)
M[, np := NP[scen]]
M[, use_text := !design_degenerate & is.finite(estimate) & abs(estimate) < 10 & is.finite(se) & se > 0 & se < 10 & (is.na(p) | is.finite(p))]
M[is.na(use_text), use_text := FALSE]
M[, use := !design_degenerate & U0 %in% TRUE]
s9_drop_const(M, "use", what = "MMM with the check")   # D-82
s11_drop_cv(M, "use", rows = rep(TRUE, nrow(M)), what = "MMM with the check")   # D-90
cat("MMM: rows", nrow(M), " use (U0)", M[, sum(use)], " use (text, no converged)", M[, sum(use_text)], " converged FALSE but finite:", M[use_text & !use, .N], "\n")
M[, z := estimate / se]; M[, rej := use & 2 * pt(-abs(z), N - np) < 0.05]
Mm <- M[, c(K6, "design_degenerate", "use", "rej", "z"), with = FALSE][, strategy := "mmm_cc"]
X <- rbind(G, Mm); rm(M, Mm, G); invisible(gc())
saveRDS(X, file.path(OUT, "per_replicate_gee_ar1_fg_mmm_cc.rds"))
tm("per-replicate table saved")

# ---- gate: per-point counts equal those behind Table 4 (exact-interval per-point files) ---------------------------
X[, arm := fifelse(label == "mmlb", "H1", "H0")]
DP <- X[arm == "H0", .(n = sum(use), kk = sum(rej)), by = .(strategy, dgp, scen, N, ntime)]
A <- rbind(fread("tables/V5_S7f_cp_S7d_dp_rates.csv"), fread("tables/V5_S7f_cp_S7e_dp_rates.csv"), fill = TRUE)
A <- unique(A[ref == "t" & strategy %in% c("gee_ar1_fg", "study_mmm", "mmm_cc")][strategy == "study_mmm", strategy := "mmm_cc"], by = c("dgp", "scen", "N", "ntime", "strategy"))
GM <- merge(DP, A[, .(strategy, dgp, scen, N, ntime, n_t = n, kk_t = kk)], by = c("strategy", "dgp", "scen", "N", "ntime"))
stopifnot(nrow(GM) == 4752L, GM[, all(n == n_t & kk == kk_t)]); cat("gate: per-point n and rejections equal Table 4 inputs at 4752 strategy-points\n")
# ---- block bootstrap over replicate indices, jointly over the cells of each model (Section 2.5) --------------------
set.seed(20261015L); B <- 2000L
X[, model := fifelse(scen == "M1c", "1c", sub("^M", "", scen))]
X[, key := paste(dgp, scen, N, ntime)]
KEYS <- unique(X[, .(key, dgp, scen, N, ntime, model)])[order(key)]
KEYS[, Nb := fifelse(N == 10, "10", fifelse(N <= 25, "15-25", ">25"))][, Tm := fifelse(ntime == 2, "2", fifelse(ntime == 3, "3", "4+"))]
mats <- function(s, arm_, what) { Y <- X[strategy == s & arm == arm_]; Mx <- matrix(0, nrow(KEYS), 1000L)
  Mx[cbind(match(Y$key, KEYS$key), Y$sim_id)] <- as.numeric(Y[[what]]); Mx }
models <- unique(KEYS$model)
Wb <- lapply(models, function(m) { Cm <- matrix(0L, 1000L, B + 1L); Cm[, 1] <- 1L; for (b in 1:B) Cm[, b + 1] <- tabulate(sample.int(1000L, 1000L, TRUE), 1000L); Cm }); names(Wb) <- models
bootdp <- function(num, den) { out <- matrix(NA_real_, nrow(KEYS), B + 1L)
  for (m in models) { r <- which(KEYS$model == m); out[r, ] <- (num[r, , drop = FALSE] %*% Wb[[m]]) / (den[r, , drop = FALSE] %*% Wb[[m]]) }; 100 * out }
# DEC-T2 (S8b, D-74): the D-11 denominator (replicates with an estimable design) comes from the design records of every
# replicate (RECS), not from the GEE rows: replicates with a degenerate design and no GEE row were counted as estimable.
RU <- unique(RECS[, c(K6, "design_degenerate"), with = FALSE], by = K6)
RU[, arm := fifelse(label == "mmlb", "H1", "H0")][, key := paste(dgp, scen, N, ntime)]
ndm <- function(arm_) { Y <- RU[arm == arm_ & key %in% KEYS$key]
  stopifnot(nrow(Y) == 1000L * nrow(KEYS), Y[, .(n = .N), by = key][, all(n == 1000L)])
  Mx <- matrix(NA_real_, nrow(KEYS), 1000L); Mx[cbind(match(Y$key, KEYS$key), Y$sim_id)] <- 1 - as.numeric(Y$design_degenerate)
  stopifnot(!anyNA(Mx)); Mx }
ND0 <- ndm("H0"); ND1 <- ndm("H1")
cat(sprintf("D-11 denominators from the design records: H0 %d, H1 %d estimable of %d; GEE-row denominators were H0 %d, H1 %d\n",
            sum(ND0), sum(ND1), 1000L * nrow(KEYS), sum(1 - mats("gee_ar1_fg", "H0", "design_degenerate")), sum(1 - mats("gee_ar1_fg", "H1", "design_degenerate"))))
RES <- list()
for (s in c("gee_ar1_fg", "mmm_cc")) RES[[s]] <- list(rate = bootdp(mats(s, "H0", "rej"), mats(s, "H0", "use")), pu = bootdp(mats(s, "H0", "use"), ND0),
                                                      upow = bootdp(mats(s, "H1", "rej"), ND1))
tm("bootstrap done")
RG <- list(all = rep(TRUE, nrow(KEYS)), T2 = KEYS$Tm == "2", T3p = KEYS$Tm != "2", N10_T3p = KEYS$N == 10 & KEYS$Tm != "2",
           N15p_T2 = KEYS$N >= 15 & KEYS$Tm == "2", N15p_T3p = KEYS$N >= 15 & KEYS$Tm != "2", small = KEYS$N <= 25 & KEYS$ntime <= 3)
for (nb in c("10", "15-25", ">25")) for (tm_ in c("2", "3", "4+")) RG[[paste0("map|", nb, "|", tm_)]] <- KEYS$Nb == nb & KEYS$Tm == tm_
for (nb in c("10", "15-25", ">25")) RG[[paste0("Nb|", nb, "|3+")]] <- KEYS$Nb == nb & KEYS$Tm != "2"
mad <- function(R, sel) colMeans(abs(R[sel, , drop = FALSE] - 5), na.rm = TRUE)
mn <- function(R, sel) colMeans(R[sel, , drop = FALSE], na.rm = TRUE)
MARGIN <- c(mad_gee_minus_mmm = 1, pu_gee_minus_mmm = 2.5, upow_mmm_minus_gee = 2)   # Section 2.5 (type I error averaged over cells, usable share, power)
cls <- function(e, s, m) { lo <- e - 2 * s; hi <- e + 2 * s
  fifelse(hi < -m, "GEE", fifelse(lo > m, "MMM", fifelse(lo >= -m & hi <= m, "similar", "about the margin"))) }
OUTB <- rbindlist(lapply(names(RG), function(r) { sel <- RG[[r]]; g <- RES$gee_ar1_fg; m <- RES$mmm_cc
  Q <- list(mad_gee = mad(g$rate, sel), mad_mmm = mad(m$rate, sel), mad_gee_minus_mmm = mad(g$rate, sel) - mad(m$rate, sel),
            rate_gee = mn(g$rate, sel), rate_mmm = mn(m$rate, sel), pu_gee = mn(g$pu, sel), pu_mmm = mn(m$pu, sel),
            pu_gee_minus_mmm = mn(g$pu, sel) - mn(m$pu, sel), upow_gee = mn(g$upow, sel), upow_mmm = mn(m$upow, sel), upow_mmm_minus_gee = mn(m$upow, sel) - mn(g$upow, sel))
  rbindlist(lapply(names(Q), function(q) { v <- Q[[q]]; data.table(region = r, n_dp = sum(sel), q = q, est = v[1], mcse = sd(v[-1])) })) }))
OUTB[q %in% names(MARGIN), margin := MARGIN[q]]
OUTB[q == "mad_gee_minus_mmm", class := cls(est, mcse, 1)]
OUTB[q %in% c("pu_gee_minus_mmm", "upow_mmm_minus_gee"), class := fifelse(abs(est) + 2 * mcse <= margin, "similar", fifelse(abs(est) - 2 * mcse > margin, "different", "about the margin"))]
fwrite(OUTB, "tables/V5_S7f_h2h_boot.csv")
print(dcast(OUTB[q %in% c("mad_gee", "mad_mmm", "mad_gee_minus_mmm")], region + n_dp ~ q, value.var = "est"), digits = 3)
print(OUTB[!is.na(class), .(region, q, est = round(est, 3), mcse = signif(mcse, 2), class)])

# ---- series for Figure 1 (panels b, d): MCSE of the two strategies by exact T and by N within T = 2 / T >= 3 ----------
SER <- list()
for (tt in sort(unique(KEYS$ntime))) SER[[paste0("T|", tt)]] <- KEYS$ntime == tt
for (nn in sort(unique(KEYS$N))) for (tb in c("2", "3+")) SER[[paste0("NT|", nn, "|", tb)]] <- KEYS$N == nn & (if (tb == "2") KEYS$ntime == 2 else KEYS$ntime >= 3)
OS <- rbindlist(lapply(names(SER), function(r) { sel <- SER[[r]]; if (!any(sel)) return(NULL)
  rbindlist(lapply(c("gee_ar1_fg", "mmm_cc"), function(s) { v <- mad(RES[[s]]$rate, sel); data.table(series = r, strategy = s, n_dp = sum(sel), mad = v[1], mad_se = sd(v[-1])) })) }))
fwrite(OS, "tables/V5_S7f_h2h_series_boot.csv"); cat("series rows", nrow(OS), "\n")
tm("done")
