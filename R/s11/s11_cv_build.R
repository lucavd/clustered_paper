# =============================================================================
# s11_cv_build.R -- S11 (D-90, D-91): variant row set "cv" = variant_nc with U0, U1 and U2 set to FALSE on the rows of
# the MMM with the convergence check (method mmm, scenarios M3-M5) whose dataset is a D-90 no-result key
# (results_v4/S11/cv_keys.rds, main-grid keys). Run with MMM_VARIANT=cv. Read-only on results_v4/variant_nc.
# Gate GC3 (chain part): only U0/U1/U2 change, only on mmm rows of M3-M5, only on keys; tables/V8_S11_cv_gate.csv.
# Outputs: results_v4/variant_cv/S4/rows_v4.rds, symlinks to the unchanged S4 inputs and study tables.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R")); source(file.path(REPO, "R/s11/s11_cv.R"))
stopifnot(identical(VARIANT, "cv"))
STUDY_TAB <- file.path(REPO, "tables"); STUDY_S4 <- file.path(REPO, "results_v4/S4")
V <- readRDS(file.path(REPO, "results_v4/variant_nc/S4/rows_v4.rds")); V0 <- copy(V)
rows <- V$method == "mmm" & V$scen %in% c("M3", "M4", "M5")
n0 <- s11_drop_cv(V, "U0", rows); n1 <- s11_drop_cv(V, "U1", rows); n2 <- s11_drop_cv(V, "U2", rows)   # D-90
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
chg <- Reduce(`|`, lapply(names(V), function(cl) !same(V[[cl]], V0[[cl]])))
cols_chg <- names(V)[sapply(names(V), function(cl) any(!same(V[[cl]], V0[[cl]])))]
kk <- s11_cv_keys(); onkey <- seq_len(nrow(V)) %in% V[kk, on = S11_K6, which = TRUE, nomatch = NULL]
CH <- V0[chg, .N, by = .(method, scen, label)]
print(CH)
GATE <- data.table(check = c("U0 switched off", "U1 switched off", "U2 switched off", "rows changed", "columns changed", "changed rows all mmm M3-M5 on keys", "rows"),
                   value = c(n0, n1, n2, sum(chg), paste(cols_chg, collapse = " "), as.integer(all(rows[chg] & onkey[chg])), nrow(V)))
fwrite(CH, file.path(STUDY_TAB, "V8_S11_cv_usable_change.csv")); fwrite(GATE, file.path(STUDY_TAB, "V8_S11_cv_gate.csv")); print(GATE)
stopifnot(all(rows[chg] & onkey[chg]), all(cols_chg %in% c("U0", "U1", "U2")))
rm(V0); invisible(gc(verbose = FALSE))
saveRDS(V, file.path(OUT_S4, "rows_v4.rds"), compress = FALSE)
for (f in setdiff(list.files(STUDY_S4), c("rows_v4.rds", "boot_idx.rds"))) {
  dst <- file.path(OUT_S4, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_S4, f), dst) }
skip <- "^(V4_perf_all\\.csv\\.gz|V4_perf_primary\\.csv|V4_S4_q1_q2\\.csv|V4_S5_.*|V4_S6.*|V5_.*|V6_.*|V7_.*|V8_.*)$"
for (f in list.files(STUDY_TAB)) {
  if (dir.exists(file.path(STUDY_TAB, f)) || grepl(skip, f)) next
  dst <- file.path(TABDIR, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_TAB, f), dst) }
cat("variant cv rows written:", nrow(V), "\n")
