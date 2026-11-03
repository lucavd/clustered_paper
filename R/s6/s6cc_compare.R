# =============================================================================
# s6cc_compare.R -- D-32: predictions P9-P12 (docs/steps/S6.md) and integrity of the variant tree.
# Reads the study row set (results_v4/S4/rows_v4.rds) and the variant one (results_v4/variant_cc/S4/rows_v4.rds)
# by explicit path. Usable and rejection exactly as in R/s5/s5_build.R (U0 and not a degenerate design;
# p < 0.05). MCSE of means over cells: block bootstrap on (model, i) with the study index matrix
# results_v4/S4/boot_idx.rds (B = 2000, seed 20261015).
# Outputs: tables/V4_S6cc_preds.csv, tables/V4_S6cc_cells.csv, tables/V4_S6cc_integrity.csv
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
stopifnot(!nzchar(VARIANT))
S6CC_VAR <- Sys.getenv("S6CC_VARIANT", "cc")   # S11 (D-90, D-91): S6CC_VARIANT=cv compares the check with the D-90 rule (variant cv) with the study fit
stopifnot(S6CC_VAR %in% c("cc", "cv"))
VS <- readRDS(file.path(REPO, "results_v4/S4/rows_v4.rds")); VC <- readRDS(file.path(REPO, sprintf("results_v4/variant_%s/S4/rows_v4.rds", S6CC_VAR)))
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
prep <- function(V) {
  X <- V[method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term), .(dgp, scen, N, ntime, label, sim_id, p, U0)]
  X[, model := as.integer(substr(scen, 2, 2))]; X[, deg := FALSE]; X[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
  X[, u := U0 & !deg]; X[, r := u & fifelse(is.na(p), FALSE, p < 0.05)]; X[, .(dgp, scen, model, N, ntime, label, sim_id, u, r)] }
S <- prep(VS); C <- prep(VC)
K <- c("dgp", "scen", "model", "N", "ntime", "label", "sim_id")
W <- merge(S, C, by = K, suffixes = c("_s", "_c")); stopifnot(nrow(W) == nrow(S), nrow(W) == nrow(C))
BIDX <- readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds")); B <- nrow(BIDX[[1]])
CNT <- lapply(BIDX, function(M) { Cm <- matrix(0, 1000L, B); for (b in seq_len(B)) Cm[, b] <- tabulate(M[b, ], 1000L); Cm })
# cell x replicate matrices for one model, then (B+1) cell rates; column 1 = all replicates once
cellrates <- function(D, ucol, rcol) {
  out <- list()
  for (m in unique(D$model)) {
    d <- D[model == m]; d[, cell := paste(dgp, scen, N, ntime, label)]; cl <- unique(d$cell)
    U <- matrix(0, length(cl), 1000L); R <- U; ix <- cbind(match(d$cell, cl), d$sim_id)
    U[ix] <- as.numeric(d[[ucol]]); R[ix] <- as.numeric(d[[rcol]])
    Wm <- cbind(1, CNT[[m]]); out[[as.character(m)]] <- list(cell = cl, rate = 100 * (R %*% Wm) / (U %*% Wm), n = rowSums(U))
  }
  list(cell = unlist(lapply(out, `[[`, "cell")), rate = do.call(rbind, lapply(out, `[[`, "rate")), n = unlist(lapply(out, `[[`, "n")))
}
judge <- function(est, mc, lo, hi) { a <- est - 2 * mc; b <- est + 2 * mc
  fifelse(a >= lo & b <= hi, "confermata", fifelse(b < lo | a > hi, "smentita", "indeterminata")) }
meanboot <- function(rt) { m <- colMeans(rt, na.rm = TRUE); c(est = m[1], mcse = sd(m[-1])) }
P <- list()
# P9: H0, T 8-10, strategy
D9 <- W[label == "manual" & ntime >= 8]; c9 <- cellrates(D9, "u_c", "r_c"); m9 <- meanboot(c9$rate)
P$P9 <- data.table(id = "P9", what = "type I, MMM with check, mean over 288 H0 cells T 8-10", n_cells = length(c9$cell), est = m9[1], mcse = m9[2], lo = 4, hi = 6)
# P10: dp clearly above 7.5 (H0, all T): strategy and study
H0 <- W[label == "manual", .(n_c = sum(u_c), k_c = sum(r_c), n_s = sum(u_s), k_s = sum(r_s)), by = .(dgp, scen, N, ntime)]
H0[, `:=`(rate_c = 100 * k_c / n_c, rate_s = 100 * k_s / n_s)]
H0[, `:=`(hi_c = rate_c - 200 * sqrt(rate_c / 100 * (1 - rate_c / 100) / n_c) > 7.5, hi_s = rate_s - 200 * sqrt(rate_s / 100 * (1 - rate_s / 100) / n_s) > 7.5)]
H0[, `:=`(lo_c = rate_c + 200 * sqrt(rate_c / 100 * (1 - rate_c / 100) / n_c) < 2.5, lo_s = rate_s + 200 * sqrt(rate_s / 100 * (1 - rate_s / 100) / n_s) < 2.5)]
P$P10 <- data.table(id = "P10", what = sprintf("H0 design points clearly above 7.5%% (of %d); study: %d", nrow(H0), H0[, sum(hi_s, na.rm = TRUE)]),
                    n_cells = nrow(H0), est = H0[, sum(hi_c, na.rm = TRUE)], mcse = 0, lo = 0, hi = 8)
# P11: power difference, H1, T 8-10
D11 <- W[label == "mmlb" & ntime >= 8]; a <- cellrates(D11, "u_c", "r_c"); b <- cellrates(D11, "u_s", "r_s"); stopifnot(identical(a$cell, b$cell))
m11 <- meanboot(a$rate - b$rate)
P$P11 <- data.table(id = "P11", what = "power difference with check - study, mean over 288 H1 cells T 8-10 (points)", n_cells = length(a$cell), est = m11[1], mcse = m11[2], lo = -2.5, hi = 1)
# P12: type I difference, H0, T 2-6
D12 <- W[label == "manual" & ntime <= 6]; a <- cellrates(D12, "u_c", "r_c"); b <- cellrates(D12, "u_s", "r_s")
m12 <- meanboot(a$rate - b$rate)
P$P12 <- data.table(id = "P12", what = "type I difference with check - study, mean over 576 H0 cells T 2-6 (points)", n_cells = length(a$cell), est = m12[1], mcse = m12[2], lo = -1, hi = 1)
PR <- rbindlist(P); PR[, outcome := fifelse(id == "P10", fifelse(est <= hi, "confermata", "smentita"), judge(est, mcse, lo, hi))]
fwrite(PR, file.path(TABDIR, "V4_S6cc_preds.csv")); print(PR)
H1 <- W[label == "mmlb", .(n_c = sum(u_c), k_c = sum(r_c), n_s = sum(u_s), k_s = sum(r_s)), by = .(dgp, scen, N, ntime)]
CE <- rbind(H0[, .(arm = "H0", dgp, scen, N, ntime, n_c, k_c, n_s, k_s, hi_c, hi_s, lo_c, lo_s)], H1[, .(arm = "H1", dgp, scen, N, ntime, n_c, k_c, n_s, k_s)], fill = TRUE)
fwrite(CE, file.path(TABDIR, "V4_S6cc_cells.csv"))
if (S6CC_VAR != "cc") { cat("S6CC_VARIANT =", S6CC_VAR, ": integrity checks of the cc variant tree not rerun (variant cv: R/s11/s11_gc3_chain.R)\n"); quit(save = "no") }   # D-90
# ---- integrity: outside MMM in M3-M5 the variant row set and performance tables equal the study ones
key <- c("dgp", "scen", "N", "ntime", "label", "method", "sim_id", "term")
oth <- function(V) V[!(method == "mmm" & scen %in% c("M3", "M4", "M5"))]
r1 <- isTRUE(all.equal(oth(VS), oth(VC), check.attributes = FALSE))
PS <- fread(file.path(REPO, "tables/V4_perf_all.csv.gz")); PC <- fread(file.path(REPO, "tables/variant_cc/V4_perf_all.csv.gz"))
# outside M3-M5 nothing may change; inside M3-M5 the other methods may change only in the intersection set
# (replicates usable for all eight methods, which depends on MMM usability)
po <- function(P) { x <- P[!(scen %in% c("M3", "M4", "M5")) | (method != "mmm" & set == "own")]; setorderv(x, intersect(names(x), c("dgp", "scen", "N", "ntime", "label", "arm", "method", "usab", "set"))); x }
r2 <- isTRUE(all.equal(po(PS), po(PC), check.attributes = FALSE))
r3 <- unname(tools::md5sum(file.path(REPO, "results_v4/S4/boot_idx.rds"))) == unname(tools::md5sum(file.path(REPO, "results_v4/variant_cc/S4/boot_idx.rds")))
IN <- data.table(check = c("row set identical outside MMM rows of M3-M5", "V4_perf_all identical outside M3-M5 and, in M3-M5, for the other methods in their own set", "boot_idx.rds identical (md5)"), pass = c(r1, r2, r3))
fwrite(IN, file.path(TABDIR, "V4_S6cc_integrity.csv")); print(IN)
