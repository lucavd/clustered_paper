# =============================================================================
# s7g_agg.R -- S7g: gates GG2-GG4 and predictions PG1-PG3 of the adaptive-quadrature run (pre-registered in docs/steps/S7g.md).
# Inputs : results_v4/S7g/rows/*.rds (this run); stored Laplace rows: results_v4/variant_cc/S4/rows_v4.rds (three study
#          mechanisms) and results_v4/S7e/rows/copula_*_manual_*.rds (Copula-AR1); block-bootstrap indices
#          results_v4/S4/boot_idx.rds (B = 2000, seed 20261015, one 2000 x 1000 matrix per model; M1c uses model 1).
# Usable : u0_flag of the study (converged, finite estimate/SE/p, |estimate| < 10, 0 < SE < 10), estimable design (D-11).
# Outputs: tables/V5_S7g_agq_dp.csv (per point x nAGQ: n, nd, k_t, k_n, rates, CP flags), tables/V5_S7g_agq_diff.csv (per point:
#          paired AGQ10 - Laplace difference, t and normal), tables/V5_S7g_agq_mech.csv (per mechanism: mean difference with
#          block-bootstrap interval, flag counts), tables/V5_S7g_agq_gates.csv, tables/V5_S7g_agq_preds.csv. Read-only on inputs.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(8L); setwd("/home/user/MMM")
t00 <- Sys.time()
cp_lo <- function(k, n) fifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
cp_up <- function(k, n) fifelse(k >= n, 1, qbeta(0.975, k + 1, pmax(n - k, 1e-12)))
NP <- c(M1 = 4L, M1c = 4L, M2 = 5L)
u0 <- function(conv, e, s, p) { w <- conv & is.finite(e) & is.finite(s) & abs(e) < 10 & s > 0 & s < 10 & is.finite(p); w[is.na(w)] <- FALSE; w }
G <- list(); gate <- function(code, ok, detail) G[[length(G) + 1]] <<- data.table(gate = code, pass = isTRUE(ok), detail = detail)
## ---- read this run ---------------------------------------------------------------------------------------------------
fl <- list.files("results_v4/S7g/rows", pattern = "\\.rds$", full.names = TRUE)
L <- lapply(fl, readRDS)
REC <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); ROW <- rbindlist(lapply(L, `[[`, "rows"), fill = TRUE); rm(L)
key <- c("dgp", "scen", "N", "ntime")
nrep <- REC[, .N, by = key]
nfit <- ROW[, .N, by = c(key, "nAGQ")]
gate("GG3", length(fl) == 1230L && nrow(nrep) == 123L && all(nrep$N == 1000L) && nrow(REC) == 123000L &&
       all(nfit[nAGQ == 1L, N] == 1000L) && all(nfit[nAGQ == 10L, N] == 1000L) && all(nfit[nAGQ == 25L, N] == 100L) && all(REC$generated),
     sprintf("%d chunk files, %d points, %d replicates, %d generation errors; fits per point: nAGQ 1 %s, 10 %s, 25 %s",
             length(fl), nrow(nrep), nrow(REC), sum(!REC$generated), paste(range(nfit[nAGQ == 1L, N]), collapse = "-"),
             paste(range(nfit[nAGQ == 10L, N]), collapse = "-"), paste(range(nfit[nAGQ == 25L, N]), collapse = "-")))
# GG3 as written requires every replicate to carry its fits; replicates whose dataset could not be generated (GLMM-AR1
# calibration, Section 2.2 of the paper) carry none. Explanation check (after the gate, not a new gate): the replicates not
# generated here are exactly those not generated in the study runs (S7d recs for the three mechanisms, S7e recs for the copula).
TK7 <- unique(REC[, ..key])[, .(chunk = 1:10), by = key]
TK7[, f := fifelse(dgp == "Copula-AR1", sprintf("results_v4/S7e/rows/copula_%s_N%d_T%d_manual_c%02d.rds", tolower(scen), N, ntime, chunk),
                   sprintf("results_v4/S7d/rows/%s_%s_N%d_T%d_manual_c%02d.rds", gsub("[^A-Za-z0-9]", "", dgp), tolower(scen), N, ntime, chunk))]
stopifnot(all(file.exists(TK7$f)))
RS <- rbindlist(lapply(seq_len(nrow(TK7)), function(j) { r <- readRDS(TK7$f[j])$recs
  data.table(TK7[j, ..key], sim_id = r$sim_id, generated = r$generated, design_degenerate = r$design_degenerate) }))
NGc <- merge(REC[, c(key, "sim_id", "generated", "design_degenerate"), with = FALSE], RS, by = c(key, "sim_id"), suffixes = c("", "_study"))
gen_same <- nrow(NGc) == 123000L && NGc[, all(generated == generated_study) && all(design_degenerate == design_degenerate_study)]
G[[1]][, detail := paste0(detail, sprintf(" | explanation: replicates not generated here %d, in the study runs %d, same set and same degenerate flags: %s",
                                         REC[generated == FALSE, .N], RS[generated == FALSE, .N], gen_same))]
G[[1]][, explained := gen_same]
ROW <- merge(ROW, REC[, c(key, "sim_id", "design_degenerate"), with = FALSE], by = c(key, "sim_id"))
ROW[, np := NP[scen]]
ROW[, use := u0(converged, estimate, se, p) & !design_degenerate]
ROW[, z := estimate / se]
ROW[, rej_t := use & 2 * pt(-abs(z), N - np) < 0.05]; ROW[, rej_n := use & 2 * pnorm(-abs(z)) < 0.05]
ROW[, singular := !is.na(D) & D < 1e-6]
## ---- GG2: Laplace refit = stored rows ---------------------------------------------------------------------------------
PTS <- unique(ROW[, ..key])
S3 <- readRDS("results_v4/variant_cc/S4/rows_v4.rds")
S3 <- S3[method == "glmer" & label == "manual" & scen %in% c("M1", "M2", "M1c")]
S3 <- merge(S3, PTS, by = key)[is.na(term) | term == "time:trt", c(key, "sim_id", "estimate", "se"), with = FALSE]
cf <- list.files("results_v4/S7e/rows", pattern = "^copula_m[12]_.*_manual_c[0-9]+\\.rds$", full.names = TRUE)
S7 <- rbindlist(lapply(cf, function(f) readRDS(f)$study[method == "glmer" & (is.na(term) | term == "time:trt")]), fill = TRUE)
S7 <- merge(S7, PTS, by = key)[, c(key, "sim_id", "estimate", "se"), with = FALSE]
OLD <- rbind(S3, S7); rm(S3); invisible(gc(verbose = FALSE))
NEW <- ROW[nAGQ == 1L & (is.na(term) | term == "time:trt"), c(key, "sim_id", "estimate", "se"), with = FALSE]
CMP <- merge(NEW, OLD, by = c(key, "sim_id"), all = TRUE, suffixes = c("_new", "_old"))
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-8)
CMP[, ok := same(estimate_new, estimate_old) & same(se_new, se_old)]
gate("GG2", all(CMP$ok) && !anyDuplicated(CMP[, c(key, "sim_id"), with = FALSE]),
     sprintf("%d replicate rows compared (stored: %d, refit: %d), %d differ beyond 1e-8; max |d est| %.2g, max |d se| %.2g",
             nrow(CMP), nrow(OLD), nrow(NEW), CMP[, sum(!ok)], CMP[, max(abs(estimate_new - estimate_old), na.rm = TRUE)],
             CMP[, max(abs(se_new - se_old), na.rm = TRUE)]))
fwrite(CMP[ok == FALSE], "tables/V5_S7g_agq_GG2_mismatch.csv")
## ---- per point x nAGQ ------------------------------------------------------------------------------------------------
DP <- ROW[, .(nd = sum(!design_degenerate), n = sum(use), k_t = sum(rej_t), k_n = sum(rej_n), n_sing = sum(singular & use),
              n_warn = sum(use & !is.na(status))), by = c(key, "nAGQ")]
DP[, `:=`(pu = 100 * n / nd, rate_t = 100 * k_t / n, rate_n = 100 * k_n / n, assess = n >= 20)]
DP[, `:=`(mc_t = sqrt(rate_t * (100 - rate_t) / n), hi_t = assess & 100 * cp_lo(k_t, n) > 7.5, lo_t = assess & 100 * cp_up(k_t, n) < 2.5,
          hi_n = assess & 100 * cp_lo(k_n, n) > 7.5)]
setorderv(DP, c(key, "nAGQ")); fwrite(DP, "tables/V5_S7g_agq_dp.csv")
## ---- GG4: nAGQ 25 vs 10 on replicates 1-100 ---------------------------------------------------------------------------
W <- dcast(ROW[nAGQ %in% c(10L, 25L) & sim_id <= 100L], dgp + scen + N + ntime + sim_id ~ nAGQ, value.var = c("use", "rej_t", "z"))
W <- W[use_10 & use_25]
agree <- W[, mean(rej_t_10 == rej_t_25)]
gate("GG4", agree >= 0.995, sprintf("%d fits usable for both; same t decision in %.4f%%; median |z25 - z10| %.2g, 99th percentile %.2g",
                                    nrow(W), 100 * agree, W[, median(abs(z_25 - z_10))], W[, quantile(abs(z_25 - z_10), 0.99)]))
## ---- paired differences AGQ10 - Laplace -------------------------------------------------------------------------------
P <- dcast(ROW[nAGQ %in% c(1L, 10L)], dgp + scen + N + ntime + sim_id ~ nAGQ, value.var = c("use", "rej_t", "rej_n"))
P[, both := use_1 & use_10]
P[, `:=`(dt = fifelse(both, as.numeric(rej_t_10) - as.numeric(rej_t_1), NA_real_), dn = fifelse(both, as.numeric(rej_n_10) - as.numeric(rej_n_1), NA_real_))]
DF <- P[, .(n_both = sum(both), diff_t = 100 * mean(dt, na.rm = TRUE), mc_t = 100 * sd(dt, na.rm = TRUE) / sqrt(sum(both)),
            diff_n = 100 * mean(dn, na.rm = TRUE), mc_n = 100 * sd(dn, na.rm = TRUE) / sqrt(sum(both)),
            disc_t = sum(dt != 0, na.rm = TRUE)), by = key]
cls <- function(e, m, a, b) fifelse(e - 2 * m >= a & e + 2 * m <= b, "confirmed", fifelse(e + 2 * m < a | e - 2 * m > b, "refuted", "indeterminate"))
DF[, pg2 := cls(diff_t, mc_t, -2.5, 2.5)]
DF <- merge(DF, DP[nAGQ == 1L, c(key, "rate_t", "hi_t", "pu"), with = FALSE], by = key)
DF <- merge(DF, DP[nAGQ == 10L, c(key, "rate_t", "hi_t", "pu"), with = FALSE], by = key, suffixes = c("_lap", "_agq"))
fwrite(DF, "tables/V5_S7g_agq_diff.csv")
## ---- per mechanism: mean difference with block bootstrap (PG1) ----------------------------------------------------------
BI <- readRDS("results_v4/S4/boot_idx.rds")
P[, model := fifelse(scen == "M2", 2L, 1L)]
boot_mean <- function(sub, col) {   # mean over points of the paired difference, resampling replicate indices jointly by model
  pts <- unique(sub[, ..key]); B <- nrow(BI[[1]]); acc <- matrix(0, B, 1); np <- nrow(pts)
  for (j in seq_len(np)) {
    x <- merge(pts[j], sub, by = key)
    x <- merge(data.table(sim_id = 1:1000), x[, c("sim_id", col, "model"), with = FALSE], by = "sim_id", all.x = TRUE)   # replicates not generated: NA
    setorder(x, sim_id); stopifnot(nrow(x) == 1000L)
    d <- x[[col]]; ok <- !is.na(d); d[!ok] <- 0; m <- if (pts$scen[j] == "M2") 2L else 1L; I <- BI[[m]]
    num <- matrix(d[I], nrow = B); den <- matrix(ok[I], nrow = B)
    acc <- acc + 100 * rowSums(num) / rowSums(den)
  }
  acc / np
}
MECH <- rbindlist(lapply(c("GenBinaryY", "GLMM-native", "GLMM-AR1", "Copula-AR1"), function(g) {
  sub <- P[dgp == g]; bt <- boot_mean(sub, "dt"); est <- DF[dgp == g, mean(diff_t)]
  bn <- boot_mean(sub, "dn")
  data.table(dgp = g, n_points = DF[dgp == g, .N], mean_diff_t = est, lo95_t = quantile(bt, 0.025, na.rm = TRUE), hi95_t = quantile(bt, 0.975, na.rm = TRUE),
             mcse_t = sd(bt, na.rm = TRUE), mean_diff_n = DF[dgp == g, mean(diff_n)], lo95_n = quantile(bn, 0.025, na.rm = TRUE), hi95_n = quantile(bn, 0.975, na.rm = TRUE),
             mean_rate_lap = DF[dgp == g, mean(rate_t_lap)], mean_rate_agq = DF[dgp == g, mean(rate_t_agq)],
             hi_lap = DF[dgp == g, sum(hi_t_lap)], hi_agq = DF[dgp == g, sum(hi_t_agq)],
             hi_both = DF[dgp == g, sum(hi_t_lap & hi_t_agq)], pu_lap = DF[dgp == g, mean(pu_lap)], pu_agq = DF[dgp == g, mean(pu_agq)],
             pg2_confirmed = DF[dgp == g, sum(pg2 == "confirmed")], pg2_indeterminate = DF[dgp == g, sum(pg2 == "indeterminate")],
             pg2_refuted = DF[dgp == g, sum(pg2 == "refuted")])
}))
MECH[, pg1 := fifelse(lo95_t >= -1 & hi95_t <= 1, "confirmed", fifelse(hi95_t < -1 | lo95_t > 1, "refuted", "indeterminate"))]
fwrite(MECH, "tables/V5_S7g_agq_mech.csv")
cop <- DF[dgp == "Copula-AR1" & hi_t_lap == TRUE]
pg3 <- if (nrow(cop)) 100 * mean(cop$hi_t_agq) else NA_real_
PR <- rbind(MECH[, .(pred = "PG1", unit = dgp, value = sprintf("%.2f [%.2f, %.2f]", mean_diff_t, lo95_t, hi95_t), outcome = pg1)],
            data.table(pred = "PG2", unit = "all points", value = sprintf("%d confirmed, %d indeterminate, %d refuted of %d", sum(DF$pg2 == "confirmed"),
                       sum(DF$pg2 == "indeterminate"), sum(DF$pg2 == "refuted"), nrow(DF)),
                       outcome = if (all(DF$pg2 == "confirmed")) "confirmed" else if (any(DF$pg2 == "refuted")) "refuted in some points" else "indeterminate in some points"),
            data.table(pred = "PG3", unit = "Copula-AR1", value = sprintf("%d of %d points clearly above with Laplace stay clearly above (%.0f%%)", sum(cop$hi_t_agq), nrow(cop), pg3),
                       outcome = if (is.na(pg3)) "not assessable" else if (pg3 >= 80) "confirmed" else "refuted"))
fwrite(PR, "tables/V5_S7g_agq_preds.csv")
GT <- rbindlist(G, fill = TRUE); fwrite(GT, "tables/V5_S7g_agq_gates.csv")
print(GT); print(MECH, digits = 3); print(PR)
print(DF[order(-abs(diff_t))][1:15, .(dgp, scen, N, ntime, n_both, diff_t, mc_t, rate_t_lap, rate_t_agq, hi_t_lap, hi_t_agq, pg2)], digits = 3)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t00, units = "mins"))))
