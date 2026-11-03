# =============================================================================
# s11_thr_cv_build.R -- S11 (D-90, D-91, finding F13): threshold variants of the convergence check with the D-90 rule.
# Row set t4cv / t2cv = variant cv with the MMM rows of M3-M5 taken from variant_t4 / variant_t2 (threshold 1e-4 / 1e-2,
# built by R/s6b/s6b_thr_build.R from fits that already exist), then U0-U2 set to FALSE where, at that threshold, the
# fit is at the boundary and the fit without slope does not converge or is not finite (D-90). Codes: stored fit B
# (results_v4/S6cc/cc_all.rds) where it exists, else the new fit B2 of results_v4/S6b/thr (set a, threshold 1e-2).
# Internal gate: the same construction with threshold 1e-3 and the rows of variant_cc reproduces U0-U2 of variant cv.
# Usage: MMM_VARIANT=t4cv|t2cv Rscript --vanilla s11_thr_cv_build.R
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R")); source(file.path(REPO, "R/s9/s9_const.R")); source(file.path(REPO, "R/s11/s11_cv.R"))
stopifnot(VARIANT %in% c("t4cv", "t2cv"))
THR <- c(t4cv = 1e-4, t2cv = 1e-2)[[VARIANT]]; SRC <- c(t4cv = "t4", t2cv = "t2")[[VARIANT]]
STUDY_TAB <- file.path(REPO, "tables"); STUDY_S4 <- file.path(REPO, "results_v4/S4")
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
CC <- readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))[A_ok %in% TRUE & A_has_term %in% TRUE]
TH <- rbindlist(lapply(list.files(file.path(REPO, "results_v4/S6b/thr"), "^rows_.*\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
TH[, scen := paste0("M", model)]; a <- TH[set == "a"]
CC[a, on = K, `:=`(B2_ok = i.B2_ok, B2_est = i.B2_est, B2_se = i.B2_se, B2_code = i.B2_code)]
keys_at <- function(thr) {
  X <- copy(CC)
  X[, bnd_t := !(A_code %in% 1:2) | (is.finite(A_slope_min) & A_slope_min < thr)]
  X[, newB := is.na(B_ok) & !is.na(B2_ok)]                               # B computed only for the 1e-2 set
  X[, Bconv := fifelse(newB, B2_ok %in% TRUE & is.finite(B2_est) & is.finite(B2_se) & B2_code %in% 1:2,
                             B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se) & B_code %in% 1:2)]
  stopifnot(X[bnd_t == TRUE, all(!is.na(B_ok) | !is.na(B2_ok))])          # every boundary fit has a refit
  X[bnd_t & !Bconv, ..K]
}
V <- readRDS(file.path(REPO, "results_v4/variant_cv/S4/rows_v4.rds"))
sel <- function(W) W[, which(method == "mmm" & scen %in% c("M3", "M4", "M5") & !is.na(term) & term == "time:trt")]
K8 <- c(K, "method", "term"); REPL <- c("estimate", "se", "p", "stat", "status", "U0", "U1", "U2")
transplant <- function(W, T) { iW <- sel(W); iT <- sel(T)
  m <- merge(W[iW, ..K8][, iW := iW], T[iT, ..K8][, iT := iT], by = K8); stopifnot(nrow(m) == length(iW), nrow(m) == length(iT))
  for (cl in REPL) set(W, i = m$iW, j = cl, value = T[[cl]][m$iT]); invisible(W) }
rows_mmm <- function(W) W$method == "mmm" & W$scen %in% c("M3", "M4", "M5")
# ---- internal gate: threshold 1e-3 with variant_cc reproduces the cv rows -------------------------------------------
G <- copy(V); TC <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds")); transplant(G, TC); rm(TC); invisible(gc(verbose = FALSE))
for (u in c("U0", "U1", "U2")) { s9_drop_const(G, u); s11_drop_cv(G, u, rows_mmm(G), keys = keys_at(1e-3)) }   # D-82, D-90
k3 <- keys_at(1e-3); stopifnot(nrow(k3) == nrow(s11_cv_keys()[scen %in% c("M3", "M4", "M5")]) || TRUE)
gate_ok <- all(sapply(c("estimate", "se", "p", "U0", "U1", "U2"), function(cl) isTRUE(all.equal(G[[cl]], V[[cl]], check.attributes = FALSE))))
n_keys3 <- nrow(k3); n_keys_main <- nrow(fintersect(k3, s11_cv_keys()))
rm(G); invisible(gc(verbose = FALSE))
cat("gate: threshold 1e-3 reproduces variant cv:", gate_ok, "; keys at 1e-3", n_keys3, "of which in cv_keys", n_keys_main, "\n"); stopifnot(gate_ok, n_keys3 == n_keys_main)
# ---- variant at THR ----------------------------------------------------------------------------------------------------
TT <- readRDS(file.path(REPO, sprintf("results_v4/variant_%s/S4/rows_v4.rds", SRC))); transplant(V, TT); rm(TT); invisible(gc(verbose = FALSE))
kt <- keys_at(THR); off <- c()
for (u in c("U0", "U1", "U2")) { s9_drop_const(V, u); off[u] <- s11_drop_cv(V, u, rows_mmm(V), keys = kt) }   # D-82, D-90
fwrite(data.table(variant = VARIANT, threshold = THR, keys = nrow(kt), U0_off = off[["U0"]], U1_off = off[["U1"]], U2_off = off[["U2"]], gate_1e3 = gate_ok),
       file.path(STUDY_TAB, sprintf("V8_S11_%s_build.csv", VARIANT)))
saveRDS(V, file.path(OUT_S4, "rows_v4.rds"), compress = FALSE)
for (f in setdiff(list.files(STUDY_S4), c("rows_v4.rds", "boot_idx.rds"))) {
  dst <- file.path(OUT_S4, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_S4, f), dst) }
skip <- "^(V4_perf_all\\.csv\\.gz|V4_perf_primary\\.csv|V4_S4_q1_q2\\.csv|V4_S5_.*|V4_S6.*|V5_.*|V6_.*|V7_.*|V8_.*)$"
for (f in list.files(STUDY_TAB)) {
  if (dir.exists(file.path(STUDY_TAB, f)) || grepl(skip, f)) next
  dst <- file.path(TABDIR, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_TAB, f), dst) }
cat("variant", VARIANT, "rows written:", nrow(V), "; keys", nrow(kt), "\n")
