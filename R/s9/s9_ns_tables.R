# =============================================================================
# s9_ns_tables.R -- S9 (D-82, D-83): the per-point tables of scenarios M6-M7 of R/s6b/s6b_ns_agg.R (lines "truths" to
# "rule safety"), recomputed from the stored rows results_v4/S6b/ns/rows_ns.rds with U0 set to FALSE on the datasets
# with a constant outcome. rows_ns.rds is read, never rewritten. Truths and denominators are read from the stored
# tables (they do not depend on usability). Fallback shares and predictions PB1-PB6 concern the marginalized model only
# (no constant-outcome usable rows) and are not rewritten.
# Outputs (overwritten): tables/V4_S6b_ns_dp.csv, V4_S6b_ns_rule_safety_dp.csv, V4_S6b_ns_rule_safety.csv.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R")); source(file.path(REPO, "R/s6b/s6b_ns_agg_lib.R"))
source(file.path(REPO, "R/s9/s9_const.R"))
source(path.expand("~/MMM/R/s11/s11_cv.R"))   # D-90 (S11, D-91): no result for the MMM with the check when its refit does not converge
suppressMessages(library(fastGHQuad))
V <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds"))
U_before <- V[, .(U0 = sum(U0 %in% TRUE)), by = .(method, label)]
s9_drop_const(V, "U0", what = "rows_ns")   # D-82
CH <- merge(U_before, V[, .(U0n = sum(U0 %in% TRUE)), by = .(method, label)], by = c("method", "label"))[, dU0 := U0 - U0n]
print(CH[dU0 > 0]); fwrite(CH, file.path(TABDIR, "V6_S9_ns_usable_change.csv"))
s11_drop_cv(V, "U0", rows = V$method == "mmm", what = "rows_ns")   # D-90
TR <- fread(file.path(TABDIR, "V4_S6b_ns_truth.csv"))
OLD <- fread(file.path(TABDIR, "V4_S6b_ns_dp.csv"))
ND <- unique(OLD[, .(dgp, scen, N, ntime, nd)]); stopifnot(nrow(ND) == 384L)
V[, deg := design_degenerate %in% TRUE]
V[TR, on = c("dgp", "scen", "ntime", "label"), `:=`(tm = i.tm, tc = i.tc)]
MARG <- c(MARGINAL_METHODS, "mmm_nocc")
V[, tok := method %in% MARG | (method %in% COND_METHODS & dgp != "GenBinaryY")]
V[, th := fifelse(method %in% MARG, tm, tc)]; V[tok == FALSE, th := NA_real_]
indicators(V)
XD <- dp_rates(V[method %in% c(PRIM, SENS, "mmm_nocc")], ND)
fwrite(XD, file.path(TABDIR, "V4_S6b_ns_dp.csv"))
RS <- rule_safety_dp(XD[method %in% c(PRIM, "mmm_nocc")])
fwrite(RS, file.path(TABDIR, "V4_S6b_ns_rule_safety_dp.csv"))
sm <- function(by) RS[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_na = sum(!assess), max_rate = max(rate[assess]),
                          n_point_above = sum(assess & rate > 7.5), pu = mean(pu), upow = mean(upow)), by = by]
RSS <- rbind(sm(c("strategy", "dgp", "scen")), sm(c("strategy", "dgp", "scen", "band")), sm("strategy")[, `:=`(dgp = "all", scen = "all")], fill = TRUE)
fwrite(RSS, file.path(TABDIR, "V4_S6b_ns_rule_safety.csv"))
cat("M6-M7 tables rewritten:", nrow(XD), nrow(RS), nrow(RSS), "\n")
