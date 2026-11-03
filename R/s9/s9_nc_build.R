# =============================================================================
# s9_nc_build.R -- S9 (D-82, D-83): variant row set "nc" = variant_kcfc with U0, U1 and U2 set to FALSE on the
# datasets with a constant outcome (results_v4/S9/const_keys.rds), for every method. Run with MMM_VARIANT=nc.
# Read-only on results_v4/variant_kcfc. Gate (GS2, chain part): only U0/U1/U2 change, only on constant keys, and
# only for methods that counted them (expected gee and gee_kc); written to tables/V6_S9_nc_gate.csv.
# Outputs: results_v4/variant_nc/S4/rows_v4.rds, symlinks to the unchanged S4 inputs and study tables.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R")); source(file.path(REPO, "R/s9/s9_const.R"))
stopifnot(identical(VARIANT, "nc"))
STUDY_TAB <- file.path(REPO, "tables"); STUDY_S4 <- file.path(REPO, "results_v4/S4")
V <- readRDS(file.path(REPO, "results_v4/variant_kcfc/S4/rows_v4.rds"))
U_before <- V[, .(U0 = sum(U0 %in% TRUE), U1 = sum(U1 %in% TRUE), U2 = sum(U2 %in% TRUE)), by = .(method, label)]
n0 <- s9_drop_const(V, "U0"); n1 <- s9_drop_const(V, "U1"); n2 <- s9_drop_const(V, "U2")   # D-82
U_after <- V[, .(U0n = sum(U0 %in% TRUE), U1n = sum(U1 %in% TRUE), U2n = sum(U2 %in% TRUE)), by = .(method, label)]
CH <- merge(U_before, U_after, by = c("method", "label"))[, `:=`(dU0 = U0 - U0n, dU1 = U1 - U1n, dU2 = U2 - U2n)]
print(CH[order(method, label)])
touched <- CH[dU0 > 0 | dU1 > 0 | dU2 > 0, unique(method)]
GATE <- data.table(check = c("U0 switched off", "U1 switched off", "U2 switched off", "methods touched", "rows", "touched only gee and gee_kc"),
                   value = c(n0, n1, n2, paste(sort(touched), collapse = " "), nrow(V), as.integer(all(touched %in% c("gee", "gee_kc")))))
fwrite(CH, file.path(STUDY_TAB, "V6_S9_nc_usable_change.csv")); fwrite(GATE, file.path(STUDY_TAB, "V6_S9_nc_gate.csv")); print(GATE)
stopifnot(all(touched %in% c("gee", "gee_kc")))
saveRDS(V, file.path(OUT_S4, "rows_v4.rds"), compress = FALSE)
for (f in setdiff(list.files(STUDY_S4), c("rows_v4.rds", "boot_idx.rds"))) {
  dst <- file.path(OUT_S4, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_S4, f), dst) }
skip <- "^(V4_perf_all\\.csv\\.gz|V4_perf_primary\\.csv|V4_S4_q1_q2\\.csv|V4_S5_.*|V4_S6.*|V5_.*|V6_.*)$"
for (f in list.files(STUDY_TAB)) {
  if (dir.exists(file.path(STUDY_TAB, f)) || grepl(skip, f)) next
  dst <- file.path(TABDIR, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_TAB, f), dst) }
cat("variant nc rows written:", nrow(V), "\n")
