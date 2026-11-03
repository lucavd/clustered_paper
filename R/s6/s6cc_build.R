# =============================================================================
# s6cc_build.R -- D-32: gates G8-G10 on the MMM reanalysis and construction of the variant row set.
# Run with MMM_VARIANT=cc (OUT_S4, TABDIR point to the variant tree; the study files are read by
# explicit path and never written).
# Inputs : results_v4/S6cc/rows/*.rds (s6_mmm_cc_run.R), results_v4/S4/rows_v4.rds (study)
# Outputs: results_v4/S6cc/cc_all.rds (per replicate), tables/V4_S6cc_gates.csv (study tables dir),
#          <variant>/S4/rows_v4.rds (study rows, MMM time:trt rows of M3-M5 replaced by the strategy,
#          U0-U2 recomputed with add_usability), symlinks to the unchanged S4 inputs and tables.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
stopifnot(identical(VARIANT, "cc"))
STUDY_TAB <- file.path(REPO, "tables"); STUDY_S4 <- file.path(REPO, "results_v4/S4")
FL <- list.files(file.path(REPO, "results_v4/S6cc/rows"), pattern = "\\.rds$", full.names = TRUE)
CC <- rbindlist(lapply(FL, readRDS), fill = TRUE)
CC[, scen := paste0("M", model)]
# A_has_term: written by the patched run script; chunks written before the patch stored a fit without the
# time:trt term as a row with estimate, SE and p all missing (data.table fill), so it is recovered from that.
if (!"A_has_term" %in% names(CC)) CC[, A_has_term := NA]
CC[A_ok == TRUE & is.na(A_has_term), A_has_term := !(is.na(A_est) & is.na(A_se) & is.na(A_p))]
saveRDS(CC, file.path(REPO, "results_v4/S6cc/cc_all.rds"))
V <- readRDS(file.path(STUDY_S4, "rows_v4.rds"))
S <- V[method == "mmm" & scen %in% c("M3", "M4", "M5")]
same <- function(a, b) (is.finite(a) & is.finite(b) & abs(a - b) < 1e-8) | (!is.finite(a) & !is.finite(b))
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
# ---- G9 completeness ---------------------------------------------------------------------------
n_files <- length(FL)
st_ok <- unique(S[!is.na(term) & term == "time:trt", ..K]); st_fail <- unique(S[is.na(term), ..K])
cc_ok <- unique(CC[A_ok == TRUE & A_has_term == TRUE, ..K]); cc_fail <- unique(CC[A_ok == FALSE, ..K])
g9 <- n_files == 17280L && nrow(fsetdiff(st_ok, cc_ok)) == 0 && nrow(fsetdiff(cc_ok, st_ok)) == 0 &&
      nrow(fsetdiff(st_fail, cc_fail)) == 0 && nrow(fsetdiff(cc_fail, st_fail)) == 0
# ---- G8 reproduction ---------------------------------------------------------------------------
M <- merge(S[!is.na(term) & term == "time:trt", c(K, "estimate", "se", "p", "status"), with = FALSE], CC[A_ok == TRUE & A_has_term == TRUE], by = K)
M[, ok := same(estimate, A_est) & same(se, A_se) & same(p, A_p)]
g8 <- M[, sum(!ok)] == 0L && nrow(M) == nrow(st_ok)
status_same <- M[, sum((is.na(status) & is.na(A_status)) | (!is.na(status) & !is.na(A_status) & status == A_status))]
# ---- G10 independent recomputation of the strategy ----------------------------------------------
pa <- strsplit(CC$A_alpha, ";", fixed = TRUE); pn <- strsplit(CC$A_alpha_names, ";", fixed = TRUE)
smin <- mapply(function(v, n) { if (length(v) == 0 || is.na(v[1])) return(NA_real_); k <- startsWith(n, "log(sigma):") & n != "log(sigma):(Intercept)"
  if (!any(k)) NA_real_ else min(exp(as.numeric(v[k]))) }, pa, pn)
sel <- CC$A_ok %in% TRUE & CC$A_has_term %in% TRUE; R <- CC[sel]; sm <- smin[sel]
R[, bnd2 := !(A_code %in% c(1L, 2L)) | (!is.na(sm) & sm < 1e-3)]
R[, useB2 := bnd2 & B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se)]
R[, `:=`(e2 = fifelse(useB2, B_est, A_est), s2 = fifelse(useB2, B_se, A_se), p2 = fifelse(useB2, B_p, A_p))]
g10 <- R[, sum(!(same(e2, cc_est) & same(s2, cc_se) & same(p2, cc_p)) | bnd2 != bnd)] == 0L
GATES <- data.table(gate = c("G8", "G9", "G10"), pass = c(g8, g9, g10),
  detail = c(sprintf("%d rows compared, %d mismatches (estimate, SE, p); status identical in %d", nrow(M), M[, sum(!ok)], status_same),
             sprintf("%d chunk files; stored time:trt keys %d, rerun %d; stored failures %d, rerun %d; fits without the term %d", n_files, nrow(st_ok), nrow(cc_ok), nrow(st_fail), nrow(cc_fail), CC[A_ok == TRUE & A_has_term == FALSE, .N]),
             sprintf("%d replicates recomputed, %d mismatches; boundary %d, fit without slope used %d", nrow(R), R[, sum(!(same(e2, cc_est) & same(s2, cc_se) & same(p2, cc_p)) | bnd2 != bnd)], R[, sum(bnd2)], R[, sum(useB2)])))
fwrite(GATES, file.path(STUDY_TAB, "V4_S6cc_gates.csv")); print(GATES)
stopifnot(all(GATES$pass))
# ---- variant row set -----------------------------------------------------------------------------
V[, rid := .I]
X <- V[method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term) & term == "time:trt", c("rid", K), with = FALSE]
X <- merge(X, CC[A_ok == TRUE & A_has_term == TRUE, c(K, "cc_est", "cc_se", "cc_p", "cc_status", "used"), with = FALSE], by = K)
stopifnot(nrow(X) == nrow(st_ok))
# status: the stored study status where the strategy keeps the study fit (it feeds U1/U2 only); the fit without
# slope gets the warnings captured during its own mm() call (warnings raised later by summary() are not captured).
X[, new_status := fifelse(used == "B", cc_status, V$status[rid])]
V[X$rid, `:=`(estimate = X$cc_est, se = X$cc_se, p = X$cc_p, stat = X$cc_est / X$cc_se, status = X$new_status)]
Y <- V[X$rid]; add_usability(Y); V[X$rid, `:=`(U0 = Y$U0, U1 = Y$U1, U2 = Y$U2)]
V[, rid := NULL]
saveRDS(V, file.path(OUT_S4, "rows_v4.rds"), compress = FALSE)
for (f in setdiff(list.files(STUDY_S4), c("rows_v4.rds", "boot_idx.rds"))) {
  dst <- file.path(OUT_S4, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_S4, f), dst) }
skip <- "^(V4_perf_all\\.csv\\.gz|V4_perf_primary\\.csv|V4_S4_q1_q2\\.csv|V4_S5_.*|V4_S6.*)$"
for (f in list.files(STUDY_TAB)) {
  if (dir.exists(file.path(STUDY_TAB, f)) || grepl(skip, f)) next
  dst <- file.path(TABDIR, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_TAB, f), dst) }
cat("variant rows written:", nrow(V), "| MMM rows replaced:", nrow(X), "| using the fit without slope:", X[, sum(used == "B")], "\n")
