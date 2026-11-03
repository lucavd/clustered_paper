# =============================================================================
# s2_ab_usable_missing.R -- S2(a) Q1 and S2(b) Q2.1-Q2.2, from the rows loaded
# by s2_load.R (results_full/, read-only).
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
A <- readRDS(file.path(OUTDIR, "rows_timetrt.rds"))
CELL <- c("dgp", "model", "N", "ntime", "label")
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink"); LRT <- "glmer_lrt"

# ---- Q2: generated replicates per cell, technique completeness ---------------
G <- A[, .(n_gen = uniqueN(sim_id)), by = CELL]
stopifnot(nrow(G) == 2880L)
ids <- unique(A[, c(CELL, "sim_id"), with = FALSE])
full <- CJ(dgp = unique(A$dgp), model = 1:5, N = N_LEVELS, ntime = T_LEVELS, label = c("mmlb", "manual"), sim_id = 1:1000)
miss <- full[!ids, on = c(CELL, "sim_id")]
fwrite(miss, file.path(TABDIR, "V4_S2_missing_ids.csv"))
fwrite(G, file.path(TABDIR, "V4_S2_replicates_by_cell.csv"))
DP <- G[, .(n_gen_min = min(n_gen), n_gen_H1 = n_gen[label == "mmlb"], n_gen_H0 = n_gen[label == "manual"]),
        by = .(dgp, model, N, ntime)]
cat(sprintf("Q2 total generated %d of %d; missing %d\n", sum(G$n_gen), 2880000L, nrow(miss)))
cat(sprintf("Q2.1 design points with < 1000 in at least one arm: %d of %d; min %d; max among those %d\n",
            sum(DP$n_gen_min < 1000), nrow(DP), min(DP$n_gen_min), max(DP[n_gen_min < 1000]$n_gen_min)))
cat("Q2.1 missing by dgp x arm:\n"); print(miss[, .N, by = .(dgp, label)])
cat("Q2.1 missing by dgp x model x ntime:\n"); print(dcast(miss[, .N, by = .(dgp, model, ntime)], dgp + model ~ ntime, value.var = "N", fill = 0))
cat(sprintf("Q2.1 cells where H0 and H1 miss the same sim_ids: %d of %d cells with any missing\n",
  { m2 <- miss[, .(ids = paste(sort(sim_id), collapse = ",")), by = .(dgp, model, N, ntime, label)]
    w <- dcast(m2, dgp + model + N + ntime ~ label, value.var = "ids"); sum(w$mmlb == w$manual, na.rm = TRUE) },
  uniqueN(miss[, .(dgp, model, N, ntime)])))
A[, tech := fifelse(method %in% c("gee_md", "gee_kc", "gee_fg", "gee_ss"), "gee_ss",
             fifelse(method == "shrink", "glmer", method))]
TC <- unique(A[, c(CELL, "sim_id", "tech"), with = FALSE])[, .(n_tech = .N), by = c(CELL, "sim_id")]
cat(sprintf("Q2.2 generated sim_ids with fewer than 5 techniques: %d (of %d)\n", sum(TC$n_tech < 5), nrow(TC)))
SS <- A[tech == "gee_ss", .(k = .N, has_null = any(method == "gee_ss")), by = c(CELL, "sim_id")]
cat(sprintf("Q2.2 gee_ss: replicates with 3 correction rows %d, with a 'fit returned NULL' row %d, other %d\n",
            sum(SS$k == 3 & !SS$has_null), sum(SS$has_null), sum(!(SS$k == 3 & !SS$has_null) & !SS$has_null)))
GL <- A[tech == "glmer", .(has_g = any(method == "glmer"), has_s = any(method == "shrink"), conv = any(method == "glmer" & converged)), by = c(CELL, "sim_id")]
cat(sprintf("Q2.2 glmer: shrink row present iff glmer converged: %s\n", all(GL$has_s == GL$conv)))

# ---- Q1: usability and warnings ------------------------------------------------
A <- A[method != "gee_ss"]           # 'fit returned NULL' rows: the 3 corrections are absent -> count as unusable
A[, U0 := fifelse(method %in% WALD,
                  converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p),
                  converged & is.finite(p))]
A[is.na(U0), U0 := FALSE]
CLS <- c(w_grad = "max\\|grad\\|", w_hess = "Hessian|unable to evaluate scaled gradient",
         w_unident = "nearly unidentifiable", w_optim = "failure to converge in|convergence code|bobyqa|Nelder_Mead",
         w_glm01 = "fitted probabilities numerically 0 or 1", w_glmnc = "algorithm did not converge",
         w_nan = "NaNs produced", w_nared = "NA/NaN replaced by maximum positive value")
W <- A[U0 == TRUE & !is.na(status)]
for (k in names(CLS)) W[, (k) := grepl(CLS[[k]], status)]
W[, w_other := !Reduce(`|`, .SD), .SDcols = names(CLS)]
W[, primary := fcase(w_grad | w_hess | w_unident | w_optim, "lme4_convergence",
                     w_glmnc, "glm_not_converged", w_glm01, "glm_fitted_0_1",
                     w_nan | w_nared, "nan_or_na_replaced", default = "other")]
cat(sprintf("Q1.2 usable rows with a warning: %d; classified: %d; 'other': %d\n",
            nrow(W), sum(W$primary != "other"), sum(W$primary == "other")))
if (any(W$w_other)) { cat("Q1.2 'other' statuses:\n"); print(W[w_other == TRUE, .N, by = .(method, status)][order(-N)]) }
fwrite(W[, c(.(n = .N), lapply(.SD, sum)), by = .(method), .SDcols = c(names(CLS), "w_other")], file.path(TABDIR, "V4_S2_warning_classes_by_method.csv"))
cat("Q1.2 primary class by method (usable rows with warning):\n"); print(dcast(W[, .N, by = .(method, primary)], method ~ primary, value.var = "N", fill = 0))

A[, warn := U0 & !is.na(status)]
A[W, on = c(CELL, "sim_id", "method"), conv_warn := i.primary %in% c("lme4_convergence", "glm_not_converged")]
A[is.na(conv_warn), conv_warn := FALSE]
U <- A[, .(n_rows = .N, n_U0 = sum(U0), n_warn = sum(warn), n_convwarn = sum(conv_warn),
           n_U1 = sum(U0 & !warn), n_U2 = sum(U0 & !conv_warn)), by = c(CELL, "method")]
# methods absent for a replicate (gee_ss NULL) add to the denominator through n_gen
U <- G[U, on = CELL]
U[, `:=`(p_U0 = n_U0 / n_gen, p_U1 = n_U1 / n_gen, p_U2 = n_U2 / n_gen,
         p_warn_among_U0 = fifelse(n_U0 > 0, n_warn / n_U0, NA_real_))]
U[, `:=`(mcse_p_U0 = mcse_p(p_U0, n_gen), mcse_p_U1 = mcse_p(p_U1, n_gen), mcse_p_U2 = mcse_p(p_U2, n_gen),
         mcse_p_warn = fifelse(n_U0 > 0, mcse_p(p_warn_among_U0, n_U0), NA_real_))]
fwrite(U, file.path(TABDIR, "V4_S2_usable_by_cell.csv"))
cat("Q1 pooled by method: share usable U0/U1/U2 and warn share among usable\n")
print(U[, .(n_gen = sum(n_gen), U0 = sum(n_U0) / sum(n_gen), U1 = sum(n_U1) / sum(n_gen), U2 = sum(n_U2) / sum(n_gen),
            warn_among_U0 = sum(n_warn) / sum(n_U0), cells_warn_gt_10pct = sum(p_warn_among_U0 > 0.10, na.rm = TRUE),
            cells_warn_gt_50pct = sum(p_warn_among_U0 > 0.50, na.rm = TRUE)), by = method][order(method)], digits = 4)
cat("Q1 warn share among usable, by method x dgp:\n")
print(dcast(U[, .(w = round(100 * sum(n_warn) / sum(n_U0), 2)), by = .(method, dgp)], method ~ dgp, value.var = "w"))

# unusable reasons (E3): non-usable rows by method and reason class
NU <- A[U0 == FALSE]
NU[, reason := fcase(!converged & grepl("timeout", status), "geesmv_timeout",
                     !converged & grepl("number of observations", status), "glmer_more_RE_than_obs",
                     !converged & grepl("singular", status), "singular_system",
                     !converged & grepl("Cgee|fitted value very close", status), "gee_fitted_value_near_1",
                     !converged & grepl("NULL", status), "fit_returned_NULL",
                     !converged, "other_error",
                     converged & (!is.finite(estimate) | !is.finite(se) | !is.finite(p)), "non_finite_estimate_se_p",
                     converged & (abs(estimate) >= 10 | se >= 10 | se <= 0), "estimate_or_se_out_of_range",
                     default = "other")]
R <- NU[, .N, by = .(method, reason)][order(method, -N)]
fwrite(R, file.path(TABDIR, "V4_S2_unusable_reasons.csv"))
cat("unusable reasons:\n"); print(dcast(R, method ~ reason, value.var = "N", fill = 0))
cat("other_error statuses (top 20):\n"); print(NU[reason == "other_error", .N, by = .(method, s = substr(status, 1, 80))][order(-N)][1:20])

# ---- Q1.1 against O_level_all (H0) and O_coverage_all (H1), set 'own' ---------
OL <- fread(file.path(TABDIR, "O_level_all.csv"))[set == "own"]
OC <- fread(file.path(TABDIR, "O_coverage_all.csv"))[set == "own"]
chk <- function(O, lab) {
  M <- merge(U[label == lab, .(dgp, model, N, ntime, method, n_U0)], O[, .(dgp, model, N, ntime, method, n)],
             by = c("dgp", "model", "N", "ntime", "method"), all = TRUE)
  M[, n_U0 := fifelse(is.na(n_U0), 0L, n_U0)]
  msub <- M[method %in% unique(O$method)]
  list(rows_in_O = nrow(O), matched = msub[!is.na(n), sum(n == n_U0)], mism = msub[!is.na(n) & n != n_U0],
       absent_nonzero = msub[is.na(n) & n_U0 > 0], absent_zero = msub[is.na(n) & n_U0 == 0][, .N])
}
for (x in list(list(OL, "manual", "O_level_all (H0)"), list(OC, "mmlb", "O_coverage_all (H1)"))) {
  r <- chk(x[[1]], x[[2]])
  cat(sprintf("Q1.1 %s: rows %d, n matched %d, mismatched %d; absent in O with n_U0 > 0: %d; absent with 0: %d\n",
              x[[3]], r$rows_in_O, r$matched, nrow(r$mism), nrow(r$absent_nonzero), r$absent_zero))
  if (nrow(r$mism)) print(head(r$mism, 20))
  if (nrow(r$absent_nonzero)) print(head(r$absent_nonzero, 20))
}
