# =============================================================================
# s6b_ns_diag.R -- S6b part B, explanation of the outcomes (post hoc, existing replicates only):
#  (1) GB3: replicates without a row for some analysis -> which analysis, and what it returned when re-run
#  (2) PB1: share of fits where the rule uses the fit without slope, by N group and series band, and why
#  (3) PB3: Fay-Graubard with t reference, null rates of all design points of the population-scale scenarios
#  (4) PB4: type I of MMM with the rule by series band and N group
# Outputs: tables/V4_S6b_ns_diag_gb3.csv, V4_S6b_ns_diag_fallback.csv, V4_S6b_ns_diag_fg.csv, V4_S6b_ns_diag_mmm.csv
# =============================================================================
source(path.expand("~/MMM/R/s6b/s6b_ns_common.R")); source(file.path(REPO, "R/s4/s4_common.R"))
OUT_NS <- file.path(REPO, "results_v4/S6b/ns")
FL <- list.files(file.path(OUT_NS, "rows"), "\\.rds$", full.names = TRUE)
L <- parallel::mclapply(FL, function(f) { x <- readRDS(f); list(r = x$rows[, .(dgp, scen, N, ntime, label, method, sim_id, term)], cc = x$cc) }, mc.cores = 16)
RW <- rbindlist(lapply(L, `[[`, "r")); CC <- rbindlist(lapply(L, `[[`, "cc"), fill = TRUE); rm(L); invisible(gc())
tech_of <- function(m) fifelse(m %in% c("glmer", "shrink"), "glmer", fifelse(m %in% c("gee_md", "gee_kc", "gee_fg", "gee_ss"), "gee_ss", m))
RW[, tech := tech_of(method)]
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
have <- RW[, .(techs = list(unique(tech))), by = K]
have[, miss := vapply(techs, function(t) paste(setdiff(TECHS_ALL, t), collapse = ","), "")]
M <- have[nzchar(miss), c(K, "miss"), with = FALSE]
cat("replicates with a missing analysis:", nrow(M), "\n"); print(M[, .N, by = miss])
# re-run the missing analyses (same seed, same generator) and record what they return
rerun <- function(r) {
  sc <- NS_SCEN[dgp == r$dgp & k == as.integer(sub("M", "", r$scen))]
  env <- load_env(r$dgp, RUNDIR_V4); gen <- make_gen(env, r$dgp, sc$slope)
  set_rep_seed(sc$k, r$sim_id); d <- as.data.frame(gen(r$N, r$ntime, r$sim_id, r$label))
  rbindlist(lapply(strsplit(r$miss, ",")[[1]], function(tc) {
    x <- tryCatch(as.data.table(env$do_analysis(d, r$N, r$ntime, r$label, 3L, tc, r$sim_id)), error = function(e) data.table(term = paste("ERROR", conditionMessage(e))))
    data.table(r[, ..K], tech = tc, n_rows = nrow(x), terms = paste(unique(x$term), collapse = ";"), has_timetrt = any(x$term %in% "time:trt"),
               events = sum(d$binY), treated_beyond_0 = sum(d$trt == 1 & d$time > 0))
  }))
}
G3 <- rbindlist(parallel::mclapply(seq_len(nrow(M)), function(j) tryCatch(rerun(M[j]), error = function(e) data.table(M[j, ..K], tech = "RERUN ERROR", terms = conditionMessage(e))), mc.cores = 16))
fwrite(G3, file.path(TABDIR, "V4_S6b_ns_diag_gb3.csv"))
print(G3[, .N, by = .(tech, has_timetrt, terms)][order(-N)])
# ---- (2) fallback share ------------------------------------------------------------------------------
U <- CC[A_ok %in% TRUE & A_has_term %in% TRUE]
U[, `:=`(band = fifelse(ntime <= 3, "T2-3", fifelse(ntime <= 6, "T4-6", "T8-10")), Ng = fifelse(N <= 20, "N10-20", fifelse(N <= 55, "N25-55", "N65-145")),
         arm = fifelse(label == "mmlb", "H1", "H0"), nonconv = !(A_code %in% c(1L, 2L)), small = !is.na(A_slope_min) & A_slope_min < 1e-3)]
FBD <- U[, .(n = .N, used_B = 100 * mean(used == "B"), nonconv = 100 * mean(nonconv), small_slope_conv = 100 * mean(small & !nonconv),
             code3 = 100 * mean(A_code %in% 3L), code4 = 100 * mean(A_code %in% 4L), code5 = 100 * mean(A_code %in% 5L),
             rej_B_H0 = 100 * mean(cc_p[used == "B"] < 0.05, na.rm = TRUE)), by = .(dgp, scen, arm, Ng, band)]
setorder(FBD, dgp, scen, arm, Ng, band); fwrite(FBD, file.path(TABDIR, "V4_S6b_ns_diag_fallback.csv"))
print(FBD[arm == "H0"], digits = 3)
# ---- (3) FG|t and MMM|t null rates, population scale ---------------------------------------------------------
RS <- fread(file.path(TABDIR, "V4_S6b_ns_rule_safety_dp.csv"))
FG <- RS[strategy %in% c("gee_fg|t", "mmm|t") & dgp == "GenBinaryY", .(scen, N, ntime, strategy, n, rate, mc, lo95 = rate - 2 * mc, hi)]
fwrite(FG, file.path(TABDIR, "V4_S6b_ns_diag_fg.csv"))
print(dcast(FG[strategy == "gee_fg|t"], scen + N ~ ntime, value.var = "rate"), digits = 3)
# ---- (4) MMM with the rule by band and N group -------------------------------------------------------------------
X <- RS[strategy %in% c("mmm|t", "mmm|n", "mmm_nocc|n", "gee_fg|t")]
X[, Ng := fifelse(N <= 20, "N10-20", fifelse(N <= 55, "N25-55", "N65-145"))]
MM <- X[, .(n_dp = .N, mean_rate = mean(rate, na.rm = TRUE), max_rate = max(rate, na.rm = TRUE), n_hi = sum(hi)), by = .(strategy, dgp, scen, band, Ng)]
setorder(MM, strategy, dgp, scen, band, Ng); fwrite(MM, file.path(TABDIR, "V4_S6b_ns_diag_mmm.csv"))
print(MM[strategy %in% c("mmm|n", "mmm|t") & dgp == "GenBinaryY" & scen == "M7"], digits = 3)
