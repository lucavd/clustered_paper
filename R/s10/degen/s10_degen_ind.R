# =============================================================================
# s10_degen_ind.R -- S10 "Degeneri", step 2a: per-replicate usable and rejection indicators of every strategy of
# Tables 3-5, rebuilt from the per-replicate rows with the definitions of the scripts that produced the tables:
#   family T5 (Table 5, 2,376 points, 19 strategies):
#     14 GEE strategies  R/s7d/s7d_agg.R (S7d rows) and R/s7e/s7e_agg.R (S7e gee rows): use = !design_degenerate &
#                        gee_ok & |est| < 10 & 0 < se < 10, minus D-82 constant keys; rej = t at N - p
#     mmm_cc             R/s7d/s7d_agg.R "study_mmm" (study rows, U0 & !degenerate, minus D-82) and s7e_agg.R mmm_cc
#     glmer, shrink      R/s7b/s7b_df.R (t at N - p) and s7e_agg.R;  glmer_lrt: p < 0.05 (chi-square reference)
#     flip               R/s7d2/s7d2_flip_rules.R "with D-57 rules": use = !design_degenerate & is.finite(p0_guard)
#                        (that script does not apply the D-82 filter; kept as in the file, reported in the join)
#   family T34 (Tables 3 and 4, main grid M1-M5 x 3 mechanisms, 8 methods): R/s5/s5_build.R on
#     results_v4/variant_nc/S4/rows_v4.rds: u0 = U0 & !degenerate (S2 table); rej = p < 0.05 (normal reference)
# For the Wald strategies also est and se (class L). Only usable rows are stored.
# Output: results_v4/S10/degen/ind/<family>_<strategy>.rds (K6 + rej + est + se), points.rds. Read-only elsewhere.
# =============================================================================
source(path.expand("~/MMM/R/s7b/s7b_common.R")); source(path.expand("~/MMM/R/s9/s9_const.R"))
suppressMessages(library(data.table)); setDTthreads(4L)
t00 <- Sys.time(); OUT <- file.path(REPO, "results_v4/S10/degen/ind"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); NPv <- NP
LOGF <- file.path(REPO, "results_v4/S10/degen/logs/ind.log")
say <- function(...) { m <- sprintf("%s %s", format(Sys.time()), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
CNT <- list()
put <- function(X, family, strategy, wald) {
  stopifnot(all(c(K6, "use", "rej") %in% names(X)))
  Y <- X[use == TRUE, c(K6, "rej", if (wald) c("est", "se")), with = FALSE]
  Y[, `:=`(N = as.integer(N), ntime = as.integer(ntime), sim_id = as.integer(sim_id))]
  stopifnot(!anyDuplicated(Y[, ..K6]))
  saveRDS(Y, file.path(OUT, sprintf("%s_%s.rds", family, strategy)))
  CNT[[length(CNT) + 1L]] <<- data.table(family = family, strategy = strategy, rows = nrow(X), usable = nrow(Y), rej = sum(Y$rej))
  say("%s %s: %d rows, %d usable, %d rejections", family, strategy, nrow(X), nrow(Y), sum(Y$rej))
}
tdf <- function(z, N, np) fifelse(is.finite(z), 2 * pt(-abs(z), N - np) < 0.05, FALSE)
## ---- 1. S7d run rows (1,896 points): 14 GEE strategies ------------------------------------------------------------
f7d <- list.files(file.path(REPO, "results_v4/S7d/rows"), pattern = "\\.rds$", full.names = TRUE); stopifnot(length(f7d) == 37920L)
L <- parallel::mclapply(f7d, function(f) { x <- readRDS(f)
  list(r = x$recs[, c(K6, "design_degenerate"), with = FALSE],
       w = x$rows[, c(K6, "cor", "est", "gee_ok", "se_unc", "se_md", "se_kc", "se_fg", "se_cr2"), with = FALSE]) }, mc.cores = 4L)
RECS <- rbindlist(lapply(L, `[[`, "r")); ROWS <- rbindlist(lapply(L, `[[`, "w")); rm(L); invisible(gc())
say("S7d: %d recs, %d rows", nrow(RECS), nrow(ROWS))
ROWS <- merge(ROWS, RECS, by = K6); ROWS[, np := NPv[scen]]
SEC <- c(unc = "se_unc", md = "se_md", kc = "se_kc", fg = "se_fg", cr2 = "se_cr2")
## ---- 2. S7e run rows (480 copula points) --------------------------------------------------------------------------------
f7e <- list.files(file.path(REPO, "results_v4/S7e/rows"), pattern = "\\.rds$", full.names = TRUE); stopifnot(length(f7e) == 9600L)
L <- parallel::mclapply(f7e, function(f) { x <- readRDS(f); x[c("recs", "study", "cc", "gee")] }, mc.cores = 4L)
ER <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); EST <- rbindlist(lapply(L, `[[`, "study"), fill = TRUE)
ECC <- rbindlist(lapply(L, `[[`, "cc"), fill = TRUE); EG <- rbindlist(lapply(L, `[[`, "gee"), fill = TRUE); rm(L); invisible(gc())
say("S7e: %d recs, %d study, %d cc, %d gee", nrow(ER), nrow(EST), nrow(ECC), nrow(EG))
EBASE <- ER[, c(K6, "design_degenerate"), with = FALSE]; EBASE[, np := NPv[scen]]
for (cn in c("exch", "ind", "ar1")) for (ty in names(SEC)) {
  if (ty == "cr2" && cn == "ar1") next
  s <- paste0("gee_", cn, "_", ty)
  X <- ROWS[cor == cn, c(K6, "est", "np", "design_degenerate", "gee_ok", SEC[[ty]]), with = FALSE]; setnames(X, SEC[[ty]], "se")
  X[, use := !design_degenerate & gee_ok & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(X, "use", what = s); X[, rej := use & tdf(est / se, N, np)]
  Y <- merge(EBASE, EG[cor == cn, c(K6, "est", "gee_ok", SEC[[ty]]), with = FALSE], by = K6, all.x = TRUE); setnames(Y, SEC[[ty]], "se")
  Y[, use := !design_degenerate & gee_ok %in% TRUE & is.finite(est) & abs(est) < 10 & is.finite(se) & se > 0 & se < 10]
  s9_drop_const(Y, "use", what = paste(s, "copula")); Y[, rej := use & tdf(est / se, N, np)]
  put(rbind(X[, c(K6, "use", "rej", "est", "se"), with = FALSE], Y[, c(K6, "use", "rej", "est", "se"), with = FALSE]), "T5", s, TRUE)
}
rm(ROWS, EG); invisible(gc())
## copula study strategies (s7e_agg.R)
S <- EST[method %in% c("mmm", "gee", "glmer", "shrink", "glmer_lrt") & term %in% "time:trt"]
Sm <- merge(S[method == "mmm"], ECC[A_ok == TRUE & A_has_term == TRUE, c(K6, "cc_est", "cc_se", "cc_p", "cc_status", "used"), with = FALSE], by = K6, all.x = TRUE)
Mcc <- copy(Sm); Mcc[!is.na(cc_est), `:=`(estimate = cc_est, se = cc_se, p = cc_p, status = fifelse(used == "B", cc_status, status))]; Mcc[, method := "mmm_cc"]
S <- rbind(S[method != "mmm"], Mcc[, names(S), with = FALSE]); S[, U0 := u0_flag(method, converged, estimate, se, p)]
COP <- list()
for (m in c("mmm_cc", "glmer", "shrink", "glmer_lrt")) {
  X <- merge(EBASE, S[method == m, c(K6, "estimate", "se", "p", "U0"), with = FALSE], by = K6, all.x = TRUE)
  X[, use := !design_degenerate & U0 %in% TRUE]; s9_drop_const(X, "use", what = paste(m, "copula"))
  if (m == "glmer_lrt") X[, rej := use & p < 0.05] else X[, rej := use & tdf(estimate / se, N, np)]
  setnames(X, "estimate", "est"); COP[[m]] <- X[, c(K6, "use", "rej", "est", "se"), with = FALSE]
}
rm(EST, ECC, S, Sm, Mcc); invisible(gc())
## ---- 3. study rows: variant_nc (M1-M5, M1c) and rows_ns (M6-M7) ------------------------------------------------------
KEEP <- c("dgp", "scen", "N", "ntime", "label", "method", "sim_id", "term", "estimate", "se", "p", "U0")
V <- readRDS(file.path(REPO, "results_v4/variant_nc/S4/rows_v4.rds")); V <- V[method %in% PRIM8 & !is.na(term), ..KEEP]; invisible(gc())
say("variant_nc rows: %d (PRIM8, term not NA); terms: %s", nrow(V), paste(unique(V$term), collapse = ","))
DEG <- fread(file.path(REPO, "tables/V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]; V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
V[, u0 := U0 & !deg]
# family T34: main grid (no M1c), the 8 methods, normal reference (stored p)
for (m in PRIM8) { X <- V[method == m & scen != "M1c"]; X[, use := u0 %in% TRUE]; X[, rej := use & fifelse(is.na(p), FALSE, p < 0.05)]
  setnames(X, "estimate", "est"); put(X, "T34", m, m != "glmer_lrt") }
B <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds")); B <- B[method %in% PRIM8 & !is.na(term), c(KEEP, "design_degenerate"), with = FALSE]; invisible(gc())
B[, deg := design_degenerate %in% TRUE]; B[, u0 := U0 & !deg]; s9_drop_const(B, "u0", what = "rows_ns")
V[, np := NPv[scen]]; B[, np := NPv[scen]]
for (m in c("mmm", "glmer", "shrink", "glmer_lrt")) {
  X <- rbind(V[method == m, .(dgp, scen, N, ntime, label, sim_id, estimate, se, p, u0, np)], B[method == m, .(dgp, scen, N, ntime, label, sim_id, estimate, se, p, u0, np)])
  X[, use := u0 %in% TRUE]
  if (m == "glmer_lrt") X[, rej := use & fifelse(is.na(p), FALSE, p < 0.05)] else X[, rej := use & tdf(estimate / se, N, np)]
  setnames(X, "estimate", "est"); sn <- if (m == "mmm") "mmm_cc" else m
  put(rbind(X[, c(K6, "use", "rej", "est", "se"), with = FALSE], COP[[sn]]), "T5", sn, m != "glmer_lrt")
}
rm(V, B, COP); invisible(gc())
## ---- 4. sign-flip test with the D-57 rules (S7d2 rows, 2,376 points) ------------------------------------------------
f2 <- list.files(file.path(REPO, "results_v4/S7d2/rows"), pattern = "\\.rds$", full.names = TRUE); stopifnot(length(f2) == 47520L)
X <- rbindlist(parallel::mclapply(f2, function(f) readRDS(f)[, .(dgp, scen, N, ntime, label, sim_id, design_degenerate, p0_guard)], mc.cores = 4L))
X[, use := !design_degenerate & is.finite(p0_guard)]; X[, rej := use & p0_guard < 0.05]
put(X, "T5", "flip", FALSE); rm(X); invisible(gc())
fwrite(rbindlist(CNT), file.path(OUT, "counts.csv"))
say("done in %.1f min", as.numeric(difftime(Sys.time(), t00, units = "mins")))
