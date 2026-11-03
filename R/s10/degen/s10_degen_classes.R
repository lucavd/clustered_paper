# =============================================================================
# s10_degen_classes.R -- S10 "Degeneri": reads the 47,520 chunk files of s10_degen_run.R and produces
#   * gate GD2: generated, design_degenerate, n_obs, sum_y equal to the stored recs for all 4,752,000 datasets; the
#     reconstructed constant-outcome datasets (generated, non-degenerate, sum_y 0 or n_obs) equal results_v4/S9/const_keys.rds;
#   * the LP vs detectseparation agreement on every dataset where detectseparation ran (all LP-separated datasets);
#   * results_v4/S10/degen/classes.rds (one row per dataset), tables/V7_S10_degen_classes.csv (counts and shares by class,
#     one by one and hierarchical, overall / mechanism / region / arm), tables/V7_S10_degen_gates.csv (GD2 part).
# Regions: "T2-3 & N<=25" (small), "N>=45", "other". Hierarchical order: Si, S without Si, V1, E1, B0, B1, none
# (Si is a subset of S, so the order of the pre-registered table, S first, would leave Si empty).
# =============================================================================
source(path.expand("~/MMM/R/s10/degen/s10_degen_common.R")); setDTthreads(4L)
TK <- dg_tasks(); stopifnot(all(file.exists(TK$out)))
t0 <- Sys.time()
L <- parallel::mclapply(TK$out, function(f) { x <- readRDS(f); list(x = x, secs = attr(x, "secs")) }, mc.cores = 4L)
X <- rbindlist(lapply(L, `[[`, "x")); SECS <- vapply(L, function(z) if (is.null(z$secs)) NA_real_ else z$secs, 0); rm(L); invisible(gc())
cat(sprintf("read %d files, %d datasets in %.1f min\n", nrow(TK), nrow(X), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
GT <- list(); gate <- function(name, ok, detail) { GT[[length(GT) + 1L]] <<- data.table(gate = "GD2", check = name, pass = isTRUE(ok), detail = detail); cat(name, ok, detail, "\n") }
eqn <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
X[, `:=`(m_gen = !eqn(generated, rec_generated), m_deg = !eqn(design_degenerate, rec_deg), m_nobs = !eqn(n_obs, rec_n_obs), m_sumy = !eqn(sum_y, rec_sum_y))]
gate("datasets", nrow(X) == 4752000L && !anyDuplicated(X[, ..K6]), sprintf("%d rows, %d unique keys", nrow(X), uniqueN(X[, ..K6])))
for (s in c("S7d", "S7e")) gate(sprintf("%s: generated, design_degenerate, n_obs, sum_y equal to the recs", s),
  X[src == s, sum(m_gen | m_deg | m_nobs | m_sumy)] == 0L,
  X[src == s, sprintf("%d datasets; mismatches generated %d, degenerate %d, n_obs %d, sum_y %d", .N, sum(m_gen), sum(m_deg), sum(m_nobs), sum(m_sumy))])
KC <- as.data.table(readRDS("/home/user/MMM/results_v4/S9/const_keys.rds"))
CK <- X[generated == TRUE & design_degenerate == FALSE & const == TRUE, ..K6]
a1 <- fsetdiff(CK, KC[, names(CK), with = FALSE]); a2 <- fsetdiff(KC[, names(CK), with = FALSE], CK)
gate("constant-outcome datasets equal the 9,711 keys of results_v4/S9/const_keys.rds", nrow(CK) == 9711L && nrow(a1) == 0L && nrow(a2) == 0L,
     sprintf("rebuilt %d, keys %d, only rebuilt %d, only keys %d; constant among degenerate designs %d", nrow(CK), nrow(KC), nrow(a1), nrow(a2),
             X[generated == TRUE & design_degenerate == TRUE & const == TRUE, .N]))
gate("non-degenerate datasets all generated", X[design_degenerate == FALSE & generated == FALSE, .N] == 0L, sprintf("%d not generated with an estimable design", X[design_degenerate == FALSE & generated == FALSE, .N]))
gate("occasion coding 0..T-1 in every generated dataset", X[generated == TRUE, all(time_ok)], sprintf("%d datasets outside", X[generated == TRUE & !time_ok, .N]))
E <- X[elig == TRUE]
LPDS <- data.table(gate = "GD3 (population part)", check = c("eligible datasets", "LP solution status not 0 (NA)", "LP separated", "detectseparation run",
  "LP separated, detectseparation not separated", "detectseparation errors", "Si (time:trt infinite)", "LP separated but time:trt finite"),
  value = c(nrow(E), E[, sum(is.na(lp))], E[, sum(S)], E[, sum(ds_run)], E[S == TRUE, sum(!(ds_outcome %in% TRUE))], E[, sum(!is.na(ds_err))], E[, sum(Si)], E[S == TRUE, sum(!Si)]))
print(LPDS)
# ---- classes.rds --------------------------------------------------------------------------------------------------
keep <- c(K6, "src", "generated", "design_degenerate", "n_obs", "sum_y", "const", "elig", "S", "Si", "V1", "E1", "B0", "B1", "lp", "ds_outcome", "ds_inf_terms",
          "n_subj", "n_var_subj", "n_ev_subj", "n_ne_subj", "ev_arm0", "ev_arm1", "nobs_arm0", "nobs_arm1", "m2_arm0", "m2_arm1",
          "occ_noev", "occ_allev", "cell_empty", "cell_noev", "cell_allev", "t_gen", "t_lp", "t_ds")
CL <- X[, ..keep]; saveRDS(CL, file.path(OUT_DG, "classes.rds"))
TIM <- data.table(what = c("chunks", "sum of chunk wall times (core-h)", "generation (core-h)", "LP (core-h)", "detectseparation (core-h)", "median chunk s", "max chunk s"),
                  value = c(length(SECS), sum(SECS, na.rm = TRUE) / 3600, sum(X$t_gen, na.rm = TRUE) / 3600, sum(X$t_lp, na.rm = TRUE) / 3600, sum(X$t_ds, na.rm = TRUE) / 3600,
                            median(SECS, na.rm = TRUE), max(SECS, na.rm = TRUE)))
fwrite(TIM, file.path(OUT_DG, "logs/timing_summary.csv")); print(TIM)
rm(X); invisible(gc())
# ---- class table ----------------------------------------------------------------------------------------------------
CL[, region := fifelse(ntime <= 3 & N <= 25, "T2-3 & N<=25", fifelse(N >= 45, "N>=45", "other"))]
CL[, arm := fifelse(label == "mmlb", "H1", "H0")]
CL[, hier := fifelse(Si %in% TRUE, "Si", fifelse(S %in% TRUE, "S not Si", fifelse(V1 %in% TRUE, "V1", fifelse(E1 %in% TRUE, "E1", fifelse(B0 %in% TRUE, "B0", fifelse(B1 %in% TRUE, "B1", "none"))))))]
summ <- function(Y) {
  e <- Y[elig == TRUE]; ne <- nrow(e)
  out <- data.table(datasets = nrow(Y), generated = sum(Y$generated), degenerate_design = sum(Y$design_degenerate),
                    constant = Y[design_degenerate == FALSE & const %in% TRUE, .N], eligible = ne)
  for (b in c("S", "Si", "V1", "E1", "B0", "B1")) { set(out, j = b, value = sum(e[[b]])); set(out, j = paste0("pct_", b), value = 100 * sum(e[[b]]) / ne) }
  for (h in c("Si", "S not Si", "V1", "E1", "B0", "B1", "none")) set(out, j = paste0("hier_", gsub(" ", "_", h)), value = sum(e$hier == h))
  out[, `:=`(any_class = ne - sum(e$hier == "none"), S_not_V1E1B0B1 = e[S & !(V1 | E1 | B0 | B1), .N], Si_not_V1E1B0B1 = e[Si & !(V1 | E1 | B0 | B1), .N],
             arm_no_subject_beyond_t0 = e[m2_arm0 == 0 | m2_arm1 == 0, .N],
             pct_occ_noev = 100 * e[, mean(occ_noev > 0)], pct_occ_allev = 100 * e[, mean(occ_allev > 0)],
             pct_cell_empty = 100 * e[, mean(cell_empty > 0)], pct_cell_noev = 100 * e[, mean(cell_noev > 0)], pct_cell_allev = 100 * e[, mean(cell_allev > 0)],
             mean_cell_noev = e[, mean(cell_noev)], S_with_no_empty_or_eventless_cell = e[S & cell_noev == 0 & cell_allev == 0 & cell_empty == 0, .N])]
  out
}
TABC <- rbind(cbind(scope = "all", mechanism = "all", region = "all", arm = "all", summ(CL)),
              cbind(scope = "mechanism", CL[, summ(.SD), by = .(mechanism = dgp)][, `:=`(region = "all", arm = "all")]),
              cbind(scope = "region", CL[, summ(.SD), by = region][, `:=`(mechanism = "all", arm = "all")]),
              cbind(scope = "arm", CL[, summ(.SD), by = arm][, `:=`(mechanism = "all", region = "all")]),
              cbind(scope = "mechanism x region", CL[, summ(.SD), by = .(mechanism = dgp, region)][, arm := "all"]),
              cbind(scope = "source", CL[, summ(.SD), by = .(mechanism = src)][, `:=`(region = "all", arm = "all")]), fill = TRUE)
setcolorder(TABC, c("scope", "mechanism", "region", "arm"))
fwrite(TABC, file.path(TB_OUT <- "/home/user/MMM/tables", "V7_S10_degen_classes.csv"))
fwrite(rbind(rbindlist(GT), LPDS[, .(gate, check, pass = NA, detail = as.character(value))]), file.path(TB_OUT, "V7_S10_degen_gates_gd2.csv"))
print(TABC[scope %in% c("all", "region"), .(scope, region, eligible, pct_S, pct_Si, pct_V1, pct_E1, pct_B0, pct_B1)])
