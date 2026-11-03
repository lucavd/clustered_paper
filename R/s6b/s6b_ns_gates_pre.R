# =============================================================================
# s6b_ns_gates_pre.R -- S6b part B, gates before the production run (D-34):
#  GB1 generator parity: (a) H1 and H0 of each new generator receive identical engine arguments except the tested
#      coefficient; (b) each new generator receives the arguments of the v4 M3 generator under the same seed except the
#      slope (sigma for GenBinaryY, D for GLMM-native), with the expected value.
#  GB2 same engine: run_chunk_ns with seed index 3 and the unmodified v4 M3 generator reproduces, replicates 1-20 of four
#      M3 cells, every study row (results_v4/variant_cc/S4/rows_v4.rds: estimate, SE, p, U0, all 9 methods) and the
#      stored rule fits (results_v4/S6cc/cc_all.rds).
#  GB4 aggregator: the S6b aggregation applied to the study M3 rows (variant cc) reproduces tables/variant_cc/
#      V4_S6_dp_rates.csv (n, nd, k, kt) and V4_S6_rule_safety_dp.csv (rate, hi, lo, pu, upow) on every M3 design point of
#      GenBinaryY and GLMM-native; the GLMM-native M3 marginal truths of tables/V4_S4_truth_by_cell.csv are reproduced.
#  GB5 bootstrap index: matrices for k = 1..7 drawn with seed 20261015 (as R/postrun/agg_v4.R); k = 1..5 identical to
#      results_v4/S4/boot_idx.rds; k = 6, 7 saved to results_v4/S6b/ns/boot_idx_67.rds.
# Output: tables/V4_S6b_ns_gates_pre.csv
# =============================================================================
source(path.expand("~/MMM/R/s6b/s6b_ns_common.R")); source(file.path(REPO, "R/s4/s4_common.R")); source(file.path(REPO, "R/s6b/s6b_ns_agg_lib.R"))
same <- function(a, b) (is.finite(a) & is.finite(b) & abs(a - b) < 1e-8) | (!is.finite(a) & !is.finite(b))
G <- list()
# ---- GB1 -------------------------------------------------------------------------------------------
capture_args <- function(dgp, gen_fun_of_env, N, ntime, label, k, i) {
  env <- load_env(dgp, RUNDIR_V4); got <- NULL
  if (dgp == "GenBinaryY") env$GenBinaryY <- function(...) { got <<- list(...); MMLB::GenBinaryY(...) }
  else { real <- env$gen_glmm_data; env$gen_glmm_data <- function(...) { got <<- list(...); real(...) } }
  gen <- gen_fun_of_env(env); set_rep_seed(k, i); gen(N, ntime, i, label); got
}
# formulas are compared as text (each call creates its own formula environment); data columns as plain vectors
strip <- function(a, drop) { a <- a[setdiff(names(a), drop)]
  a <- lapply(a, function(x) if (inherits(x, "formula")) paste(deparse(x), collapse = "") else if (is.data.frame(x)) lapply(x, as.vector) else x); a }
res1 <- list()
for (s in seq_len(nrow(NS_SCEN))) for (N in c(10L, 55L)) for (Tt in c(3L, 10L)) for (i in c(1L, 7L)) {
  sc <- NS_SCEN[s]; bname <- if (sc$dgp == "GenBinaryY") "beta" else "beta_c"; sname <- if (sc$dgp == "GenBinaryY") "sigma" else "D"
  newg <- function(env) make_gen(env, sc$dgp, sc$slope); oldg <- function(env) env$generate_data_model3
  h1 <- capture_args(sc$dgp, newg, N, Tt, "mmlb", sc$k, i); h0 <- capture_args(sc$dgp, newg, N, Tt, "manual", sc$k, i)
  o1 <- capture_args(sc$dgp, oldg, N, Tt, "mmlb", sc$k, i)
  a_ok <- identical(strip(h1, bname), strip(h0, bname)) && identical(h1[[bname]][-4], h0[[bname]][-4]) && h0[[bname]][4] == 0 && h1[[bname]][4] > 0
  b_ok <- identical(strip(h1, sname), strip(o1, sname))
  exp_s <- if (sc$dgp == "GenBinaryY") c(1, sc$slope) else matrix(c(1, 0.01, 0.01, sc$slope), 2, byrow = TRUE)
  s_ok <- isTRUE(all.equal(unname(h1[[sname]]), unname(exp_s))) && isTRUE(all.equal(unname(h0[[sname]]), unname(exp_s)))
  res1[[length(res1) + 1]] <- data.table(dgp = sc$dgp, k = sc$k, N = N, ntime = Tt, i = i, a_ok = a_ok, b_ok = b_ok, s_ok = s_ok)
}
R1 <- rbindlist(res1); print(R1[, .(all_a = all(a_ok), all_b = all(b_ok), all_s = all(s_ok)), by = .(dgp, k)])
G$GB1 <- data.table(gate = "GB1", pass = R1[, all(a_ok & b_ok & s_ok)], detail = sprintf("%d comparisons (4 scenarios x N 10, 55 x T 3, 10 x 2 seeds): H1 vs H0 ok %d, new vs M3 ok %d, slope value ok %d", nrow(R1), R1[, sum(a_ok)], R1[, sum(b_ok)], R1[, sum(s_ok)]))
# ---- GB2 -------------------------------------------------------------------------------------------
VC <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds")); CC <- readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))
VS <- readRDS(file.path(REPO, "results_v4/S4/rows_v4.rds"))
cells <- data.table(dgp = c("GenBinaryY", "GenBinaryY", "GLMM-native", "GLMM-native"), N = c(25L, 25L, 55L, 55L), ntime = c(4L, 4L, 8L, 8L), label = c("manual", "mmlb", "manual", "mmlb"))
res2 <- list()
for (r in seq_len(nrow(cells))) {
  t <- cbind(cells[r], data.table(k = 3L, scen = "M3", slope = NA_real_, chunk = 1L))
  x <- run_chunk_ns(t, k_seed = 3L, gen_override = function(env) env$generate_data_model3, ids = 1:20, out_file = NULL)
  R <- conv_rows(copy(x$rows))
  key <- c("method", "sim_id")
  # study rows: MMM compared with the study fit (VS), other methods with the variant (identical outside MMM)
  ST <- VS[dgp == t$dgp & scen == "M3" & N == t$N & ntime == t$ntime & label == t$label & sim_id <= 20]
  M <- merge(R[, .(method, sim_id, term, estimate, se, p, U0)], ST[, .(method, sim_id, term, estimate, se, p, U0)], by = key, suffixes = c("", ".s"), all = TRUE)
  ok_rows <- M[, sum(!(same(estimate, estimate.s) & same(se, se.s) & same(p, p.s) & U0 %in% U0.s))]
  C2 <- merge(x$cc, CC[dgp == t$dgp & scen == "M3" & N == t$N & ntime == t$ntime & label == t$label & sim_id <= 20], by = "sim_id", suffixes = c("", ".s"))
  ok_cc <- C2[, sum(!(same(cc_est, cc_est.s) & same(cc_se, cc_se.s) & same(cc_p, cc_p.s) & used == used.s))]
  mA <- merge(x$cc[, .(sim_id, A_est)], R[method == "mmm" & !is.na(term), .(sim_id, estimate)], by = "sim_id")
  ok_A <- mA[, sum(!same(A_est, estimate))]
  res2[[r]] <- data.table(t[, .(dgp, N, ntime, label)], rows_compared = nrow(M), rows_mismatch = ok_rows, rows_missing = M[is.na(estimate) & !is.na(estimate.s) | is.na(U0) | is.na(U0.s), .N],
                          cc_compared = nrow(C2), cc_mismatch = ok_cc, A_vs_mmm_mismatch = ok_A)
}
R2 <- rbindlist(res2); print(R2)
G$GB2 <- data.table(gate = "GB2", pass = R2[, all(rows_mismatch == 0 & rows_missing == 0 & cc_mismatch == 0 & A_vs_mmm_mismatch == 0 & cc_compared == 20)],
  detail = sprintf("%d study rows and %d rule fits compared on 4 M3 cells x 20 replicates: %d and %d mismatches; fit A vs study MMM row %d mismatches",
                   R2[, sum(rows_compared)], R2[, sum(cc_compared)], R2[, sum(rows_mismatch + rows_missing)], R2[, sum(cc_mismatch)], R2[, sum(A_vs_mmm_mismatch)]))
# ---- GB4 -------------------------------------------------------------------------------------------
TT <- fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))
tr_chk <- TT[dgp == "GLMM-native" & scen == "M3"]
tr_chk[, mine := mapply(function(lab, T) native_truth3(if (lab == "mmlb") c(-2.32, 0.116, 0.232, 0.232) else c(-2.32, 0.116, 0.232, 0), matrix(c(1, 0.01, 0.01, 0.001), 2), T), label, ntime)]
tr_ok <- tr_chk[, max(abs(mine - truth_marginal))] < 1e-10
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[model == 3, .(N, ntime, sim_id)]
V <- VC[scen == "M3" & dgp %in% c("GenBinaryY", "GLMM-native") & !is.na(term) & method %in% c(PRIM, SENS), .(dgp, scen, N, ntime, label, method, sim_id, estimate, se, p, U0)]
V[, deg := FALSE]; V[DEG, on = c("N", "ntime", "sim_id"), deg := TRUE]
V[TT, on = c("dgp", "scen", "ntime", "label"), `:=`(tm = i.truth_marginal, tc = i.truth_conditional)]
V[, tok := method %in% MARGINAL_METHODS | (method %in% COND_METHODS & dgp != "GenBinaryY")]
V[, th := fifelse(method %in% MARGINAL_METHODS, tm, tc)]; V[tok == FALSE, th := NA_real_]
indicators(V)
ND <- CJ(dgp = c("GenBinaryY", "GLMM-native"), scen = "M3", N = as.integer(N_LEVELS), ntime = as.integer(T_LEVELS))
DC <- DEG[, .(ndeg = .N), by = .(N, ntime)]; ND[DC, on = c("N", "ntime"), ndeg := i.ndeg]; ND[is.na(ndeg), ndeg := 0L]; ND[, nd := 1000L - ndeg][, ndeg := NULL]
X <- dp_rates(V[method %in% PRIM], ND)
REF <- fread(file.path(REPO, "tables/variant_cc/V4_S6_dp_rates.csv"))[scen == "M3" & dgp %in% c("GenBinaryY", "GLMM-native")]
m4 <- merge(X, REF, by = c("dgp", "scen", "N", "ntime", "arm", "method"), suffixes = c("", ".r"))
d4 <- m4[, sum(n != n.r | nd != nd.r | k != k.r | !((is.na(kt) & is.na(kt.r)) | kt %in% kt.r))]
RS <- rule_safety_dp(X); RSR <- fread(file.path(REPO, "tables/variant_cc/V4_S6_rule_safety_dp.csv"))[scen == "M3" & dgp %in% c("GenBinaryY", "GLMM-native")]
m5 <- merge(RS, RSR, by = c("dgp", "scen", "N", "ntime", "strategy"), suffixes = c("", ".r"))
d5 <- m5[, sum(!(same(rate, rate.r) & hi == hi.r & lo == lo.r & same(pu, pu.r) & same(upow, upow.r)))]
G$GB4 <- data.table(gate = "GB4", pass = nrow(m4) == nrow(REF) && nrow(REF) == 2L * 96L * 2L * 8L && d4 == 0L && nrow(m5) == nrow(RSR) && d5 == 0L && tr_ok,
  detail = sprintf("dp rates: %d of %d rows matched, %d mismatches; rule safety: %d of %d rows, %d mismatches; GLMM-native M3 truths max |diff| %.2g",
                   nrow(m4), nrow(REF), d4, nrow(m5), nrow(RSR), d5, tr_chk[, max(abs(mine - truth_marginal))]))
# ---- GB5 -------------------------------------------------------------------------------------------
set.seed(20261015L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
BI <- lapply(1:7, function(m) matrix(sample.int(1000L, 2000L * 1000L, replace = TRUE), nrow = 2000L))
BS <- readRDS(file.path(REPO, "results_v4/S4/boot_idx.rds"))
gb5 <- all(vapply(1:5, function(m) identical(BI[[m]], BS[[m]]), logical(1)))
dir.create(OUT_NS, recursive = TRUE, showWarnings = FALSE); saveRDS(list(`6` = BI[[6]], `7` = BI[[7]]), file.path(OUT_NS, "boot_idx_67.rds"))
G$GB5 <- data.table(gate = "GB5", pass = gb5, detail = sprintf("k = 1..5 identical to results_v4/S4/boot_idx.rds: %s; k = 6, 7 saved", gb5))
GT <- rbindlist(G); fwrite(GT, file.path(TABDIR, "V4_S6b_ns_gates_pre.csv")); print(GT)
