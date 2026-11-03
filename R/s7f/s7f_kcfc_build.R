# =============================================================================
# s7f_kcfc_build.R -- D-60: variant row set "kcfc" = variant_cc with the gee_kc rows carrying the
# Kauermann-Carroll correction of fewclust 0.2.0 (S7d run, exchangeable GEE) instead of geesmv's.
# Run with MMM_VARIANT=kcfc. Read-only on results_v4/variant_cc and results_v4/S7d.
# Gate GF3(b): only gee_kc rows change; their estimate equals the study's where both exist.
# Outputs: results_v4/variant_kcfc/S4/rows_v4.rds, symlinks to the unchanged S4 inputs and tables,
#          tables/V5_S7f_kcfc_gate.csv (study tables dir).
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
stopifnot(identical(VARIANT, "kcfc"))
STUDY_TAB <- file.path(REPO, "tables"); STUDY_S4 <- file.path(REPO, "results_v4/S4")
CC_S4 <- file.path(REPO, "results_v4/variant_cc/S4")
V <- readRDS(file.path(CC_S4, "rows_v4.rds"))
V0 <- copy(V)
FL <- list.files(file.path(REPO, "results_v4/S7d/rows"), pattern = "\\.rds$", full.names = TRUE)
rd <- function(f) { x <- readRDS(f)$rows; x[cor == "exch", .(dgp, scen, N, ntime, label, sim_id, gee_ok, est, se_kc, msg_kc)] }
F <- rbindlist(lapply(FL, rd))
F <- F[scen %in% c("M1", "M2", "M3", "M4", "M5", "M1c")]
F[scen != "M1c", scen := toupper(scen)]
K <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
stopifnot(F[, .(cnt = .N), by = K][, all(cnt == 1L)])
V[, rid := .I]
X <- V[method == "gee_kc" & !is.na(term) & term == "time:trt", c("rid", K, "estimate", "se", "p", "stat", "converged", "status"), with = FALSE]
G <- V[method == "gee" & !is.na(term) & term == "time:trt", c(K, "status", "converged"), with = FALSE]
setnames(G, c("status", "converged"), c("gee_status", "gee_conv"))
X <- merge(X, F, by = K, all.x = TRUE)
X <- merge(X, G, by = K, all.x = TRUE)
n_kc <- nrow(X); n_nomatch <- X[is.na(gee_ok), .N]
# study p of gee_kc is the normal two-sided p of estimate/se (checked on the rows where it is finite)
chk <- X[is.finite(se) & is.finite(p)]
p_rule_ok <- chk[, max(abs(p - 2 * pnorm(-abs(estimate / se))))] < 1e-10
both <- X[is.finite(estimate) & is.finite(est)]
est_diff <- both[, max(abs(estimate - est))]
# new values
X[, ok_fc := gee_ok %in% TRUE & is.finite(se_kc) & se_kc > 0]
X[, `:=`(new_se = fifelse(ok_fc, se_kc, NA_real_))]
X[, `:=`(new_stat = estimate / new_se)]
X[, `:=`(new_p = 2 * pnorm(-abs(new_stat)), new_conv = ok_fc,
         new_status = fifelse(ok_fc, gee_status, fifelse(is.na(msg_kc), "fewclust: no KC standard error", paste0("fewclust: ", msg_kc))))]
V[X$rid, `:=`(se = X$new_se, stat = X$new_stat, p = X$new_p, converged = X$new_conv, status = X$new_status)]
Y <- V[X$rid]; add_usability(Y); V[X$rid, `:=`(U0 = Y$U0, U1 = Y$U1, U2 = Y$U2)]
V[, rid := NULL]
# GF3(b): nothing else changed
other <- V0$method != "gee_kc" | is.na(V0$term) | V0$term != "time:trt"
same_col <- function(a, b) { (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b) }
cols <- names(V0)
n_other_diff <- sum(!Reduce(`&`, lapply(cols, function(cc) same_col(V0[[cc]][other], V[[cc]][other]))))
unchanged_kc_cols <- setdiff(cols, c("se", "stat", "p", "converged", "status", "U0", "U1", "U2"))
kcr <- which(!other)
n_kc_other_col_diff <- sum(!Reduce(`&`, lapply(unchanged_kc_cols, function(cc) same_col(V0[[cc]][kcr], V[[cc]][kcr]))))
GATE <- data.table(
  check = c("rows identical outside gee_kc time:trt", "gee_kc columns other than se/stat/p/converged/status/U identical",
            "gee_kc estimate = S7d exchangeable estimate (max abs diff)", "study gee_kc p = normal p of estimate/se",
            "gee_kc keys without an S7d row", "gee_kc rows", "usable (U0) study KC", "usable (U0) fewclust KC",
            "S7d exch rows read"),
  value = c(n_other_diff, n_kc_other_col_diff, est_diff, as.integer(p_rule_ok), n_nomatch, n_kc,
            V0[method == "gee_kc", sum(U0)], V[method == "gee_kc", sum(U0)], nrow(F)))
GATE[, pass := c(value[1] == 0, value[2] == 0, value[3] < 1e-8, value[4] == 1, value[5] == 0, TRUE, TRUE, TRUE, TRUE)]
fwrite(GATE, file.path(STUDY_TAB, "V5_S7f_kcfc_gate.csv")); print(GATE)
stopifnot(all(GATE$pass))
saveRDS(V, file.path(OUT_S4, "rows_v4.rds"), compress = FALSE)
for (f in setdiff(list.files(STUDY_S4), c("rows_v4.rds", "boot_idx.rds"))) {
  dst <- file.path(OUT_S4, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_S4, f), dst) }
skip <- "^(V4_perf_all\\.csv\\.gz|V4_perf_primary\\.csv|V4_S4_q1_q2\\.csv|V4_S5_.*|V4_S6.*|V5_.*)$"
for (f in list.files(STUDY_TAB)) {
  if (dir.exists(file.path(STUDY_TAB, f)) || grepl(skip, f)) next
  dst <- file.path(TABDIR, f); if (!file.exists(dst)) file.symlink(file.path(STUDY_TAB, f), dst) }
cat("variant kcfc rows written:", nrow(V), "| gee_kc rows replaced:", n_kc, "\n")
