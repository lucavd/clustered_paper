# =============================================================================
# s6b_thr_compare.R -- S6b part A: predictions PA1-PA6 and gate GA5 (docs/steps/S6b.md). Variants t4 (threshold 1e-4)
# and t2 (1e-2) against cc (1e-3, the primary analysis). Definitions as R/s6/s6cc_compare.R (usable = U0 and not a
# degenerate design; rejection p < 0.05; block bootstrap on (model, i) with results_v4/S4/boot_idx.rds, B = 2000).
# Outputs: tables/V4_S6b_thr_preds.csv, V4_S6b_thr_cells.csv, V4_S6b_thr_rule.csv, V4_S6b_thr_claims_diff.csv,
#          V4_S6b_thr_integrity.csv
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
stopifnot(!nzchar(VARIANT))
rows <- function(tag) readRDS(file.path(REPO, sprintf("results_v4/variant_%s/S4/rows_v4.rds", tag)))
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
prep <- function(V) {
  X <- V[method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term), .(dgp, scen, N, ntime, label, sim_id, p, U0)]
  X[, model := as.integer(substr(scen, 2, 2))]; X[, deg := FALSE]; X[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
  X[, u := U0 & !deg]; X[, r := u & fifelse(is.na(p), FALSE, p < 0.05)]; X[, .(dgp, scen, model, N, ntime, label, sim_id, u, r)] }
K <- c("dgp", "scen", "model", "N", "ntime", "label", "sim_id")
VS <- readRDS(file.path(REPO, "results_v4/S4/rows_v4.rds"))
C0 <- prep(rows("cc"))
BIDX <- readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds")); B <- nrow(BIDX[[1]])
CNT <- lapply(BIDX, function(M) { Cm <- matrix(0, 1000L, B); for (b in seq_len(B)) Cm[, b] <- tabulate(M[b, ], 1000L); Cm })
cellrates <- function(D, ucol, rcol) {
  out <- list()
  for (m in unique(D$model)) {
    d <- D[model == m]; d[, cell := paste(dgp, scen, N, ntime, label)]; cl <- unique(d$cell)
    U <- matrix(0, length(cl), 1000L); R <- U; ix <- cbind(match(d$cell, cl), d$sim_id)
    U[ix] <- as.numeric(d[[ucol]]); R[ix] <- as.numeric(d[[rcol]])
    Wm <- cbind(1, CNT[[m]]); out[[as.character(m)]] <- list(cell = cl, rate = 100 * (R %*% Wm) / (U %*% Wm))
  }
  list(cell = unlist(lapply(out, `[[`, "cell")), rate = do.call(rbind, lapply(out, `[[`, "rate")))
}
judge <- function(est, mc, lo, hi) { a <- est - 2 * mc; b <- est + 2 * mc
  fifelse(a >= lo & b <= hi, "confermata", fifelse(b < lo | a > hi, "smentita", "indeterminata")) }
meanboot <- function(rt) { m <- colMeans(rt, na.rm = TRUE); c(est = m[1], mcse = sd(m[-1])) }
oth <- function(V) V[!(method == "mmm" & scen %in% c("M3", "M4", "M5"))]
EXP <- list(t4 = c(PA1 = 0.1, PA2 = 0.4, PA3 = 0.1), t2 = c(PA1 = 0, PA2 = -0.1, PA3 = -0.1))
PR <- list(); CE <- list(); IN <- list(); RU <- list(); CL <- list()
CLc <- fread(file.path(REPO, "tables/variant_cc/V4_S5_claims.csv"))
RSc <- fread(file.path(REPO, "tables/variant_cc/V4_S6_rule_safety.csv"))
for (tg in c("t4", "t2")) {
  VT <- rows(tg)
  IN[[tg]] <- data.table(variant = tg, check = "row set identical to the study outside MMM rows of M3-M5", pass = isTRUE(all.equal(oth(VS), oth(VT), check.attributes = FALSE)))
  C1 <- prep(VT); rm(VT); gc(verbose = FALSE)
  W <- merge(C0, C1, by = K, suffixes = c("_c", "_v")); stopifnot(nrow(W) == nrow(C0), nrow(W) == nrow(C1))
  dd <- function(D) { a <- cellrates(D, "u_v", "r_v"); b <- cellrates(D, "u_c", "r_c"); stopifnot(identical(a$cell, b$cell)); list(cell = a$cell, d = a$rate - b$rate, v = a$rate, c = b$rate) }
  x1 <- dd(W[label == "manual"]); x2 <- dd(W[label == "manual" & ntime >= 8]); x3 <- dd(W[label == "mmlb"])
  m1 <- meanboot(x1$d); m2 <- meanboot(x2$d); m3 <- meanboot(x3$d)
  P <- data.table(variant = tg, id = c("PA1", "PA2", "PA3"),
    what = c("type I difference vs cc, mean over H0 cells M3-M5 (points)", "same, H0 cells with T 8-10", "power difference vs cc, mean over H1 cells M3-M5 (points)"),
    n_cells = c(length(x1$cell), length(x2$cell), length(x3$cell)), expected = EXP[[tg]], est = c(m1[1], m2[1], m3[1]), mcse = c(m1[2], m2[2], m3[2]), lo = -1, hi = 1)
  P[, outcome := judge(est, mcse, lo, hi)]
  # PA4: per-cell paired differences, H0
  cd <- data.table(variant = tg, cell = x1$cell, d = x1$d[, 1], mc = apply(x1$d[, -1, drop = FALSE], 1, sd, na.rm = TRUE), rate_v = x1$v[, 1], rate_c = x1$c[, 1])
  cd[, outcome := judge(d, mc, -2.5, 2.5)]; CE[[tg]] <- cd
  P <- rbind(P, data.table(variant = tg, id = "PA4", what = "H0 cells M3-M5 with paired difference clearly beyond +/- 2.5 points (smentita)", n_cells = nrow(cd),
                           expected = 0, est = cd[, sum(outcome == "smentita")], mcse = 0, lo = 0, hi = 8,
                           outcome = fifelse(cd[, sum(outcome == "smentita")] <= 8, "confermata", "smentita")))
  RS <- fread(file.path(REPO, sprintf("tables/variant_%s/V4_S6_rule_safety.csv", tg)))
  r1 <- RS[is.na(N) & band == "all" & strategy %in% c("mmm|t", "mmm|n", "gee_fg|t")][, variant := tg]
  RU[[tg]] <- r1
  hi_t <- r1[strategy == "mmm|t", n_hi]
  P <- rbind(P, data.table(variant = tg, id = "PA5", what = sprintf("design points with MMM|t clearly above 7.5%% (of %d); MMM|n: %d", r1[strategy == "mmm|t", n_dp], r1[strategy == "mmm|n", n_hi]),
                           n_cells = r1[strategy == "mmm|t", n_dp], expected = 0, est = hi_t, mcse = 0, lo = 0, hi = 0, outcome = fifelse(hi_t == 0, "confermata", "smentita")))
  CLv <- fread(file.path(REPO, sprintf("tables/variant_%s/V4_S5_claims.csv", tg)))
  M <- merge(CLc[, .(name, class_cc = class, est_cc = est, lo_cc = lo, hi_cc = hi)], CLv[, .(name, class_v = class, est_v = est, lo_v = lo, hi_v = hi)], by = "name")
  M[, changed := class_cc != class_v]; CL[[tg]] <- M[changed == TRUE][, variant := tg]
  P <- rbind(P, data.table(variant = tg, id = "PA6", what = sprintf("sub-claims changing class vs cc (of %d)", nrow(M)), n_cells = nrow(M),
                           expected = 0, est = M[, sum(changed)], mcse = 0, lo = 0, hi = 0, outcome = fifelse(M[, sum(changed)] == 0, "confermata", "smentita")))
  PR[[tg]] <- P
}
RU[["cc"]] <- RSc[is.na(N) & band == "all" & strategy %in% c("mmm|t", "mmm|n", "gee_fg|t")][, variant := "cc"]
PRt <- rbindlist(PR); print(PRt[, .(variant, id, expected, est = round(est, 3), mcse = signif(mcse, 2), outcome)])
fwrite(PRt, file.path(TABDIR, "V4_S6b_thr_preds.csv")); fwrite(rbindlist(CE), file.path(TABDIR, "V4_S6b_thr_cells.csv"))
fwrite(rbindlist(RU, use.names = TRUE), file.path(TABDIR, "V4_S6b_thr_rule.csv")); fwrite(rbindlist(CL, fill = TRUE), file.path(TABDIR, "V4_S6b_thr_claims_diff.csv"))
INt <- rbindlist(IN); fwrite(INt, file.path(TABDIR, "V4_S6b_thr_integrity.csv")); print(INt)
print(rbindlist(RU, use.names = TRUE)[, .(variant, strategy, n_hi, n_lo, pu, upow, max_rate)])
print(rbindlist(CL, fill = TRUE))
print(rbindlist(CE)[, .N, by = .(variant, outcome)])
