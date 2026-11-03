# =============================================================================
# s6b_thr_build.R -- S6b part A: gates GA3-GA4 and the variant row sets with the "at the boundary" threshold
# of the D-32 rule moved to 1e-4 (variant t4) and 1e-2 (variant t2). Same construction as R/s6/s6cc_build.R.
# Inputs : results_v4/S6cc/cc_all.rds, results_v4/S6b/thr/rows_*.rds (s6b_thr_run.R), results_v4/S4/rows_v4.rds,
#          results_v4/variant_cc/S4/rows_v4.rds (GA4 only)
# Outputs: tables/V4_S6b_thr_gates.csv, tables/V4_S6b_thr_inventory.csv,
#          results_v4/variant_<tag>/S4/rows_v4.rds + symlinks to the unchanged S4 inputs and study tables
# Usage  : Rscript --vanilla s6b_thr_build.R      (no MMM_VARIANT; both variants are written)
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
stopifnot(!nzchar(VARIANT))
STUDY_TAB <- file.path(REPO, "tables"); STUDY_S4 <- file.path(REPO, "results_v4/S4")
same <- function(a, b) (is.finite(a) & is.finite(b) & abs(a - b) < 1e-8) | (!is.finite(a) & !is.finite(b))
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
CC <- readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))
TH <- rbindlist(lapply(list.files(file.path(REPO, "results_v4/S6b/thr"), "^rows_.*\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
TH[, scen := paste0("M", model)]
U <- CC[A_ok %in% TRUE & A_has_term %in% TRUE]
# ---- GA3 same engine -----------------------------------------------------------------------------
a <- TH[set == "a"]; b <- TH[set == "b"]
ga3a <- a[, sum(!(same(A2_est, A_est) & same(A2_se, A_se) & same(A2_p, A_p)))]
ga3b <- b[, sum(!(B2_ok %in% TRUE & same(B2_est, B_est) & same(B2_se, B_se)))]
# ---- inventory ------------------------------------------------------------------------------------
conv <- U$A_code %in% c(1L, 2L); sm <- U$A_slope_min
n_t4 <- sum(conv & !is.na(sm) & sm >= 1e-4 & sm < 1e-3)
n_t2 <- sum(conv & !is.na(sm) & sm >= 1e-3 & sm < 1e-2)
keyA <- U[conv & !is.na(sm) & sm >= 1e-3 & sm < 1e-2, ..K]
inv_ok <- n_t4 == 21879L && n_t2 == 19962L && nrow(a) == n_t2 && nrow(fsetdiff(keyA, a[, ..K])) == 0L && nrow(b) == 2000L
# ---- strategy with threshold thr ------------------------------------------------------------------
U[a, on = K, `:=`(B2_ok = i.B2_ok, B2_est = i.B2_est, B2_se = i.B2_se, B2_p = i.B2_p, B2_status = i.B2_status)]
strategy <- function(thr) {
  X <- U[, c(K, "A_est", "A_se", "A_p", "A_status", "A_code", "A_slope_min", "B_ok", "B_est", "B_se", "B_p", "B_status", "B2_ok", "B2_est", "B2_se", "B2_p", "B2_status"), with = FALSE]
  X[, bnd_t := !(A_code %in% c(1L, 2L)) | (!is.na(A_slope_min) & A_slope_min < thr)]
  # fit without slope: stored one where it exists (boundary at 1e-3), the new one for the set computed here
  X[, newB := !is.na(B2_ok)]
  X[, `:=`(Bo = fifelse(newB, B2_ok %in% TRUE, B_ok %in% TRUE), Be = fifelse(newB, B2_est, B_est), Bs = fifelse(newB, B2_se, B_se),
           Bp = fifelse(newB, B2_p, B_p), Bst = fifelse(newB, B2_status, B_status))]
  X[, useB := bnd_t & Bo & is.finite(Be) & is.finite(Bs)]
  X[, `:=`(est = fifelse(useB, Be, A_est), se = fifelse(useB, Bs, A_se), p = fifelse(useB, Bp, A_p), used = fifelse(useB, "B", "A"),
           st = fifelse(useB, Bst, A_status))]
  X
}
S3 <- strategy(1e-3)
M4 <- merge(S3[, c(K, "est", "se", "p", "used"), with = FALSE], U[, c(K, "cc_est", "cc_se", "cc_p", "used"), with = FALSE], by = K, suffixes = c("", ".u"))
stopifnot(nrow(M4) == nrow(U))
ga4a <- M4[, sum(!(same(est, cc_est) & same(se, cc_se) & same(p, cc_p)) | used != used.u)]
V <- readRDS(file.path(STUDY_S4, "rows_v4.rds"))
build_rows <- function(X) {
  W <- copy(V); W[, rid := .I]
  R <- W[method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term) & term == "time:trt", c("rid", K), with = FALSE]
  R <- merge(R, X[, c(K, "est", "se", "p", "st", "used"), with = FALSE], by = K)
  stopifnot(nrow(R) == nrow(X))
  R[, new_status := fifelse(used == "B", st, W$status[rid])]
  W[R$rid, `:=`(estimate = R$est, se = R$se, p = R$p, stat = R$est / R$se, status = R$new_status)]
  Y <- W[R$rid]; add_usability(Y); W[R$rid, `:=`(U0 = Y$U0, U1 = Y$U1, U2 = Y$U2)]
  W[, rid := NULL]; W
}
W3 <- build_rows(S3); VC <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))
ga4b <- isTRUE(all.equal(W3, VC, check.attributes = FALSE)); rm(W3, VC); gc(verbose = FALSE)
INV <- data.table(item = c("usable MMM fits of M3-M5", "converged, 1e-4 <= min exp(a_k) < 1e-3 (back to the study fit with 1e-4)",
                           "converged, 1e-3 <= min exp(a_k) < 1e-2 (fit without slope needed with 1e-2)", "  of which the fit without slope is usable",
                           "fit without slope used, threshold 1e-4", "fit without slope used, threshold 1e-3", "fit without slope used, threshold 1e-2"),
                  n = c(nrow(U), n_t4, n_t2, a[, sum(B2_ok %in% TRUE & is.finite(B2_est) & is.finite(B2_se))],
                        strategy(1e-4)[, sum(useB)], S3[, sum(useB)], strategy(1e-2)[, sum(useB)]))
fwrite(INV, file.path(STUDY_TAB, "V4_S6b_thr_inventory.csv")); print(INV)
GATES <- data.table(gate = c("GA3a", "GA3b", "GA4"), pass = c(ga3a == 0L, ga3b == 0L, ga4a == 0L && ga4b && inv_ok),
  detail = c(sprintf("%d study fits redone, %d mismatches (estimate, SE, p)", nrow(a), ga3a),
             sprintf("%d stored fits without slope redone, %d mismatches (estimate, SE)", nrow(b), ga3b),
             sprintf("threshold 1e-3: %d fits, %d mismatches with cc_est/cc_se/cc_p/used; row set identical to variant_cc: %s; inventory 21879/19962 found: %s (%d/%d)",
                     nrow(S3), ga4a, ga4b, inv_ok, n_t4, n_t2)))
fwrite(GATES, file.path(STUDY_TAB, "V4_S6b_thr_gates.csv")); print(GATES)
stopifnot(all(GATES$pass))
skip <- "^(V4_perf_all\\.csv\\.gz|V4_perf_primary\\.csv|V4_S4_q1_q2\\.csv|V4_S5_.*|V4_S6.*)$"
for (tg in c("t4", "t2")) {
  thr <- c(t4 = 1e-4, t2 = 1e-2)[[tg]]
  OD <- file.path(REPO, "results_v4", paste0("variant_", tg), "S4"); TD <- file.path(REPO, "tables", paste0("variant_", tg))
  dir.create(OD, recursive = TRUE, showWarnings = FALSE); dir.create(TD, recursive = TRUE, showWarnings = FALSE)
  W <- build_rows(strategy(thr))
  saveRDS(W, file.path(OD, "rows_v4.rds"), compress = FALSE); rm(W); gc(verbose = FALSE)
  for (f in setdiff(list.files(STUDY_S4), c("rows_v4.rds", "boot_idx.rds"))) { dst <- file.path(OD, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_S4, f), dst) }
  for (f in list.files(STUDY_TAB)) { if (dir.exists(file.path(STUDY_TAB, f)) || grepl(skip, f)) next
    dst <- file.path(TD, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_TAB, f), dst) }
  cat("variant", tg, "written\n")
}
