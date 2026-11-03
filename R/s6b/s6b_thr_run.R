# =============================================================================
# s6b_thr_run.R -- S6b part A: fits needed to move the "at the boundary" threshold of the D-32 rule to 0.01.
# Set (a): fits of results_v4/S6cc/cc_all.rds that converge (nlm code 1 or 2) and have 0.001 <= min exp(a_k) < 0.01,
#          i.e. not at the boundary with 0.001 but at the boundary with 0.01: refit A (study fit, gate GA3a) and B
#          (lv.formula = ~1). Set (b): 2,000 fits already at the boundary with a stored fit B, sampled with seed
#          20261019: refit B (gate GA3b). Data regenerated exactly as R/s6/s6_mmm_cc_run.R (same seed, generators,
#          formulas, q = 20, warnings captured).
# Usage: Rscript --vanilla s6b_thr_run.R <dgp> <cores>
# Output: results_v4/S6b/thr/rows_<dgp>.rds
# =============================================================================
args <- commandArgs(TRUE); DGP <- args[1]; CORES <- as.integer(args[2])
source(path.expand("~/MMM/R/s6/s6_common.R"))
suppressMessages({ library(MMLB); library(dplyr) })
data.table::setDTthreads(1L)
RD <- file.path(REPO, "R/v4/rundir"); sp <- function(f) source(file.path(RD, f))
if (DGP == "GenBinaryY") {
  for (f in c(sprintf("generation_model%d.R", 1:5), "generate_data.R", "make_all_data.R")) sp(f)
} else if (DGP == "GLMM-native") {
  sp("gen_glmm_data.R"); for (f in c(sprintf("manual/generation_model%d.R", 1:5), "manual/generate_data.R", "manual/make_all_data.R")) sp(f)
} else {
  sp("gen_ar1_data.R"); for (f in c("ar1/generation_models_ar1.R", "manual/generate_data.R", "ar1/make_all_data.R")) sp(f)
}
FORM <- list(`3` = list(m = binY ~ time * trt, lv = ~time, t = ~1),
             `4` = list(m = binY ~ time * trt, lv = ~time, t = ~trt),
             `5` = list(m = binY ~ time * trt + dose, lv = ~time + dose, t = ~trt + steroids))
fitw <- function(d, F, lv) {
  warns <- character(0)
  x <- withCallingHandlers(tryCatch(mm(F$m, lv.formula = lv, t.formula = F$t, data = d, id = id, q = 20), error = function(e) NULL),
                           warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(x)) return(list(ok = FALSE))
  s <- summary(x)$mean.table; k <- which(rownames(s) == "time:trt")
  if (!length(k)) return(list(ok = TRUE, has_term = FALSE))
  list(ok = TRUE, has_term = TRUE, est = unname(s[k, 1]), se = unname(s[k, 2]), p = unname(s[k, 4]), code = unname(x$control["convergence_code"]),
       status = if (length(warns)) paste(warns, collapse = " | ") else NA_character_)
}
CC <- readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))
U <- CC[A_ok %in% TRUE & A_has_term %in% TRUE]
setA <- U[dgp == DGP & A_code %in% c(1L, 2L) & !is.na(A_slope_min) & A_slope_min >= 1e-3 & A_slope_min < 1e-2]
set.seed(20261019L); allB <- U[bnd %in% TRUE & B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se)]
pick <- sort(sample.int(nrow(allB), 2000L)); setB <- allB[pick][dgp == DGP]   # same 2,000 drawn in every DGP call
J <- rbind(setA[, .(set = "a", model, N, ntime, label, sim_id, A_est, A_se, A_p, B_est, B_se)],
           setB[, .(set = "b", model, N, ntime, label, sim_id, A_est, A_se, A_p, B_est, B_se)])
cat(DGP, ": set a", nrow(setA), "| set b", nrow(setB), "\n")
J[, grp := .GRP, by = .(model, N, ntime, label)]
G <- split(seq_len(nrow(J)), J$grp)
run_grp <- function(ix) {
  out <- vector("list", length(ix))
  for (q in seq_along(ix)) {
    r <- J[ix[q]]; F <- FORM[[as.character(r$model)]]
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 1000L * r$model + r$sim_id)
    d <- tryCatch(as.data.frame(make_all_data(r$N, r$ntime, r$sim_id, r$model, r$label)), error = function(e) NULL)
    if (is.null(d)) { out[[q]] <- data.table(row = ix[q], gen_ok = FALSE); next }
    A <- if (r$set == "a") fitw(d, F, F$lv) else list(ok = NA)
    B <- fitw(d, F, ~1); okB <- isTRUE(B$ok) && isTRUE(B$has_term)
    out[[q]] <- data.table(row = ix[q], gen_ok = TRUE,
      A2_est = if (isTRUE(A$has_term)) A$est else NA_real_, A2_se = if (isTRUE(A$has_term)) A$se else NA_real_, A2_p = if (isTRUE(A$has_term)) A$p else NA_real_,
      B2_ok = okB, B2_est = if (okB) B$est else NA_real_, B2_se = if (okB) B$se else NA_real_, B2_p = if (okB) B$p else NA_real_,
      B2_code = if (okB) B$code else NA_integer_, B2_status = if (okB) B$status else NA_character_)
  }
  rbindlist(out, fill = TRUE)
}
t0 <- Sys.time()
R <- rbindlist(parallel::mclapply(G, function(ix) tryCatch(run_grp(ix), error = function(e) data.table(row = ix, gen_ok = NA, err = conditionMessage(e))),
                                  mc.cores = CORES, mc.preschedule = FALSE), fill = TRUE)
X <- cbind(J[R$row], R[, -"row"]); X[, dgp := DGP]
dir.create(file.path(REPO, "results_v4/S6b/thr"), recursive = TRUE, showWarnings = FALSE)
saveRDS(X, file.path(REPO, sprintf("results_v4/S6b/thr/rows_%s.rds", gsub("[^A-Za-z0-9]", "", DGP))))
cat(DGP, "done:", nrow(X), "rows,", X[is.na(gen_ok) | gen_ok == FALSE, .N], "not regenerated;", format(difftime(Sys.time(), t0, units = "mins")), "\n")
