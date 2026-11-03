# =============================================================================
# s6b_ns_agg.R -- S6b part B: gate GB3, measures and predictions PB1-PB6 on the new scenarios M6-M7 (D-34).
# Input : results_v4/S6b/ns/rows/*.rds (s6b_ns_run.R), results_v4/S6b/ns/boot_idx_67.rds (s6b_ns_gates_pre.R)
# Output: results_v4/S6b/ns/rows_ns.rds (rows_v4 layout; "mmm" = with the D-32 rule, "mmm_nocc" = study fit)
#         tables/V4_S6b_ns_gates_post.csv, V4_S6b_ns_truth.csv, V4_S6b_ns_dp.csv, V4_S6b_ns_rule_safety_dp.csv,
#         V4_S6b_ns_rule_safety.csv, V4_S6b_ns_fallback.csv, V4_S6b_ns_preds.csv
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R")); source(file.path(REPO, "R/s6b/s6b_ns_agg_lib.R"))
suppressMessages(library(fastGHQuad))
OUT_NS <- file.path(REPO, "results_v4/S6b/ns")
N_LEVELS <- c(10, 15, 20, seq(25, 145, 10)); T_LEVELS <- c(2, 3, 4, 6, 8, 10)
FL <- list.files(file.path(OUT_NS, "rows"), "\\.rds$", full.names = TRUE)
TEST <- Sys.getenv("S6B_TEST") == "1"   # code check on a subset of chunks; writes only to a temporary directory
if (TEST) { FL <- FL[unique(round(seq(1, length(FL), length.out = 300)))]; TABDIR <- tempfile("s6b_test_"); dir.create(TABDIR); OUT_SAVE <- TABDIR } else OUT_SAVE <- OUT_NS
L <- parallel::mclapply(FL, readRDS, mc.cores = 24)
RW <- rbindlist(lapply(L, `[[`, "rows"), fill = TRUE); CC <- rbindlist(lapply(L, `[[`, "cc"), fill = TRUE)
RC <- rbindlist(lapply(seq_along(L), function(j) { x <- L[[j]]$recs; x[, file := basename(FL[j])]; x }), fill = TRUE); rm(L); invisible(gc(verbose = FALSE))
RC[, c("dgp_t", "m", "N", "ntime", "label", "chunk") := tstrsplit(sub("\\.rds$", "", file), "_")]
RC[, `:=`(dgp = fifelse(dgp_t == "GenBinaryY", "GenBinaryY", "GLMM-native"), scen = toupper(m), N = as.integer(sub("N", "", N)), ntime = as.integer(sub("T", "", ntime)))]
CELLK <- c("dgp", "scen", "N", "ntime", "label")
# ---- GB3 completeness ------------------------------------------------------------------------------
tech_of <- function(m) fifelse(m %in% c("glmer", "shrink"), "glmer", fifelse(m %in% c("gee_md", "gee_kc", "gee_fg", "gee_ss"), "gee_ss", m))
RW[, tech := tech_of(method)]
cov_t <- RW[, .(nt = uniqueN(tech[tech != "generation"]), gen_fail = any(tech == "generation")), by = c(CELLK, "sim_id")]
n_rep <- RC[, .(nrep = .N), by = CELLK]
gb3 <- length(FL) == 7680L && nrow(n_rep) == 768L && all(n_rep$nrep == 1000L) &&
       nrow(cov_t) == 768000L && cov_t[, all(nt == 5L | gen_fail)] && nrow(CC) == 768000L
GP <- data.table(gate = "GB3", pass = gb3, detail = sprintf("%d chunk files; %d cells; %d replicates with rows, %d with all five analyses, %d generation failures; %d rule records",
                 length(FL), nrow(n_rep), nrow(cov_t), cov_t[nt == 5L, .N], cov_t[gen_fail == TRUE, .N], nrow(CC)))
fwrite(GP, file.path(TABDIR, "V4_S6b_ns_gates_post.csv")); print(GP)
# ---- rows in the rows_v4 layout, MMM with the rule -----------------------------------------------------
RW[, tech := NULL]
V <- conv_rows(RW[, .(dgp, scen, N = as.integer(N), ntime = as.integer(ntime), label, method, sim_id = as.integer(sim_id), term, estimate, se, stat, p, converged, status, design_degenerate)])
rm(RW); invisible(gc(verbose = FALSE))
NO <- V[method == "mmm"][, method := "mmm_nocc"]
CCu <- CC[A_ok %in% TRUE & A_has_term %in% TRUE]
stopifnot(nrow(CCu) == V[method == "mmm" & !is.na(term) & term == "time:trt", .N])
V[, rid := .I]
R <- merge(V[method == "mmm" & !is.na(term) & term == "time:trt", c("rid", CELLK, "sim_id", "status"), with = FALSE],
           CCu[, c(CELLK, "sim_id", "cc_est", "cc_se", "cc_p", "cc_status", "used", "A_est"), with = FALSE], by = c(CELLK, "sim_id"))
stopifnot(nrow(R) == nrow(CCu))
chkA <- V[R$rid, sum(!(abs(estimate - R$A_est) < 1e-8 | (!is.finite(estimate) & !is.finite(R$A_est))))]
R[, new_status := fifelse(used == "B", cc_status, status)]
V[R$rid, `:=`(estimate = R$cc_est, se = R$cc_se, p = R$cc_p, stat = R$cc_est / R$cc_se, status = R$new_status)]
Y <- V[R$rid]; add_usability(Y); V[R$rid, `:=`(U0 = Y$U0, U1 = Y$U1, U2 = Y$U2)]; V[, rid := NULL]
V <- rbindlist(list(V, NO), use.names = TRUE); rm(NO, Y); invisible(gc(verbose = FALSE))
cat("fit A of the rule vs MMM row of the study fit, mismatches:", chkA, "\n")
saveRDS(V, file.path(OUT_SAVE, "rows_ns.rds"), compress = FALSE)
# ---- truths ------------------------------------------------------------------------------------------
SL <- c(M6 = 0.01, M7 = 0.04)
TR <- rbind(CJ(dgp = "GenBinaryY", scen = c("M6", "M7"), ntime = as.integer(T_LEVELS), label = c("mmlb", "manual"))[, `:=`(tm = fifelse(label == "mmlb", 0.2, 0), tc = NA_real_)],
            CJ(dgp = "GLMM-native", scen = c("M6", "M7"), ntime = as.integer(T_LEVELS), label = c("mmlb", "manual"))[, tc := fifelse(label == "mmlb", 0.232, 0)], fill = TRUE)
TR[dgp == "GLMM-native", tm := mapply(function(s, lab, T) native_truth3(if (lab == "mmlb") c(-2.32, 0.116, 0.232, 0.232) else c(-2.32, 0.116, 0.232, 0),
                                                                        matrix(c(1, 0.01, 0.01, SL[[s]]), 2), T), scen, label, ntime)]
fwrite(TR, file.path(TABDIR, "V4_S6b_ns_truth.csv"))
V[, deg := design_degenerate %in% TRUE]
V[TR, on = c("dgp", "scen", "ntime", "label"), `:=`(tm = i.tm, tc = i.tc)]
MARG <- c(MARGINAL_METHODS, "mmm_nocc")
V[, tok := method %in% MARG | (method %in% COND_METHODS & dgp != "GenBinaryY")]
V[, th := fifelse(method %in% MARG, tm, tc)]; V[tok == FALSE, th := NA_real_]
indicators(V)
ND <- RC[label == "manual", .(nd = sum(!(design_degenerate %in% TRUE))), by = .(dgp, scen, N, ntime)]
ND2 <- RC[label == "mmlb", .(nd2 = sum(!(design_degenerate %in% TRUE))), by = .(dgp, scen, N, ntime)]
stopifnot(identical(merge(ND, ND2)[, nd], merge(ND, ND2)[, nd2]))
XD <- dp_rates(V[method %in% c(PRIM, SENS, "mmm_nocc")], ND)
fwrite(XD, file.path(TABDIR, "V4_S6b_ns_dp.csv"))
RS <- rule_safety_dp(XD[method %in% c(PRIM, "mmm_nocc")])
fwrite(RS, file.path(TABDIR, "V4_S6b_ns_rule_safety_dp.csv"))
sm <- function(by) RS[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_na = sum(!assess), max_rate = max(rate[assess]),
                          n_point_above = sum(assess & rate > 7.5), pu = mean(pu), upow = mean(upow)), by = by]
RSS <- rbind(sm(c("strategy", "dgp", "scen")), sm(c("strategy", "dgp", "scen", "band")), sm("strategy")[, `:=`(dgp = "all", scen = "all")], fill = TRUE)
fwrite(RSS, file.path(TABDIR, "V4_S6b_ns_rule_safety.csv"))
print(RSS[is.na(band) & strategy %in% c("mmm|t", "mmm|n", "mmm_nocc|t", "gee_fg|t", "gee_md|t", "glmer_lrt|n")][order(strategy, dgp, scen)], digits = 3)
# ---- fallback share and behaviour ----------------------------------------------------------------------
CCu[, band := fifelse(ntime <= 3, "T2-3", fifelse(ntime <= 6, "T4-6", "T8-10"))]
CCu[, arm := fifelse(label == "mmlb", "H1", "H0")]
FB <- rbind(CCu[, .(n = .N, at_boundary = mean(bnd), used_B = mean(used == "B"), nonconv = mean(!(A_code %in% c(1L, 2L))), small_slope = mean(A_slope_min < 1e-3, na.rm = TRUE),
                    rej_B = mean(cc_p[used == "B"] < 0.05, na.rm = TRUE), rej_A_kept = mean(cc_p[used == "A"] < 0.05, na.rm = TRUE)), by = .(dgp, scen, arm)][, band := "all"],
            CCu[, .(n = .N, at_boundary = mean(bnd), used_B = mean(used == "B"), nonconv = mean(!(A_code %in% c(1L, 2L))), small_slope = mean(A_slope_min < 1e-3, na.rm = TRUE),
                    rej_B = mean(cc_p[used == "B"] < 0.05, na.rm = TRUE), rej_A_kept = mean(cc_p[used == "A"] < 0.05, na.rm = TRUE)), by = .(dgp, scen, arm, band)])
fwrite(FB, file.path(TABDIR, "V4_S6b_ns_fallback.csv")); print(FB[band == "all"], digits = 3)
# ---- predictions ----------------------------------------------------------------------------------------
BI <- readRDS(file.path(OUT_NS, "boot_idx_67.rds")); CNT <- lapply(BI, boot_counts)
D <- dcast(V[method %in% c("mmm", "mmm_nocc"), .(dgp, scen, N, ntime, label, sim_id, method, u0, r)], dgp + scen + N + ntime + label + sim_id ~ method, value.var = c("u0", "r"))
D[, model := as.integer(sub("M", "", scen))]
P <- list()
EXP1 <- c(`GenBinaryY M6` = 15, `GenBinaryY M7` = 5, `GLMM-native M6` = 35, `GLMM-native M7` = 15)
for (g in c("GenBinaryY", "GLMM-native")) for (s in c("M6", "M7")) {
  key <- paste(g, s)
  f <- FB[dgp == g & scen == s & arm == "H0" & band == "all", 100 * used_B]
  P[[length(P) + 1]] <- data.table(id = "PB1", scenario = key, what = "share of MMM fits using the fit without slope, H0 (%)", n_cells = 96L, expected = EXP1[[key]], est = f, mcse = NA_real_, lo = 0, hi = 40)
  d0 <- D[dgp == g & scen == s & label == "manual"]; d1 <- D[dgp == g & scen == s & label == "mmlb"]
  a <- cellrates(d0, "u0_mmm", "r_mmm", CNT); b <- cellrates(d0, "u0_mmm_nocc", "r_mmm_nocc", CNT)
  a1 <- cellrates(d1, "u0_mmm", "r_mmm", CNT); b1 <- cellrates(d1, "u0_mmm_nocc", "r_mmm_nocc", CNT)
  m4 <- meanboot(a$rate); m5 <- meanboot(a$rate - b$rate); m6 <- meanboot(a1$rate - b1$rate)
  P[[length(P) + 1]] <- data.table(id = c("PB4", "PB5", "PB6"), scenario = key,
    what = c("type I of MMM with the rule, mean over the 96 H0 cells (%)", "type I difference with rule - without, mean over H0 cells (points)", "power difference with rule - without, mean over H1 cells (points)"),
    n_cells = 96L, expected = c(5.3, if (g == "GenBinaryY") -0.3 else -0.5, 0), est = c(m4[1], m5[1], m6[1]), mcse = c(m4[2], m5[2], m6[2]), lo = c(4, -1, -1), hi = c(6, 1, 1))
}
hi_m <- RS[strategy == "mmm|t", sum(hi)]; hi_f <- RS[strategy == "gee_fg|t", sum(hi)]
P[[length(P) + 1]] <- data.table(id = c("PB2", "PB3"), scenario = "all four", what = c(sprintf("design points with MMM|t clearly above 7.5%% (of %d)", RS[strategy == "mmm|t", .N]),
                                 sprintf("design points with FG|t clearly above 7.5%% (of %d)", RS[strategy == "gee_fg|t", .N])), n_cells = 384L, expected = 0, est = c(hi_m, hi_f), mcse = 0, lo = 0, hi = 0)
PR <- rbindlist(P, use.names = TRUE)
PR[, outcome := fifelse(id %in% c("PB1", "PB2", "PB3"), fifelse(est >= lo & est <= hi, "confermata", "smentita"), judge(est, mcse, lo, hi))]
fwrite(PR, file.path(TABDIR, "V4_S6b_ns_preds.csv")); print(PR[, .(id, scenario, expected, est = round(est, 3), mcse = signif(mcse, 2), outcome)])
