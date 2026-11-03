# =============================================================================
# s6_mmm_cc_run.R -- reanalysis of MMM in M3-M5 with the convergence check of D-32.
# Same datasets as the study: for each (dgp, model, N, ntime, label) and replicate i = 1..1000 the data
# are regenerated with set.seed(20260805 + 1000 * model + i), L'Ecuyer-CMRG, by the v4 generators
# (as R/v4/rundir/run_full.R). Fit A = study fit (formulas of R/v4/rundir/analysis_model3-5.R, q = 20,
# defaults; warnings captured as in do_analysis.R). "At the boundary": nlm code not in {1, 2} or some
# latent-SD slope exp(a_k) < 0.001. Then fit B = same model with lv.formula = ~1. Strategy (cc):
# B if A is at the boundary and B gives finite estimate and SE, otherwise A.
# Usage: Rscript s6_mmm_cc_run.R <dgp> <cores>      (resumable: existing chunk files are skipped)
# Output: results_v4/S6cc/rows/<dgp>_m<model>_N<N>_T<T>_<label>_c<chunk>.rds (100 replicates each)
# =============================================================================
args <- commandArgs(TRUE); DGP <- args[1]; CORES <- as.integer(args[2])
source(path.expand("~/MMM/R/s6/s6_common.R"))
suppressMessages({ library(MMLB); library(dplyr) })
data.table::setDTthreads(1L)   # one thread per worker: the default OpenMP pool (half the cores) oversubscribed the host
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
OUTD <- file.path(REPO, "results_v4/S6cc/rows"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
N_LEVELS <- c(10, 15, 20, seq(25, 145, 10)); T_LEVELS <- c(2, 3, 4, 6, 8, 10); CHUNK <- 100L
TASKS <- CJ(model = 3:5, N = as.integer(N_LEVELS), ntime = as.integer(T_LEVELS), label = c("mmlb", "manual"), chunk = 1:10)
TASKS[, file := file.path(OUTD, sprintf("%s_m%d_N%d_T%d_%s_c%02d.rds", gsub("[^A-Za-z0-9]", "", DGP), model, N, ntime, label, chunk))]
TASKS <- TASKS[!file.exists(file)][order(-N * ntime)]
if (nzchar(Sys.getenv("S6CC_MAXTASKS"))) TASKS <- head(TASKS, as.integer(Sys.getenv("S6CC_MAXTASKS")))   # smoke test only
cat(DGP, "tasks to run:", nrow(TASKS), "\n")
fitw <- function(d, F, lv) {
  warns <- character(0)
  x <- withCallingHandlers(tryCatch(mm(F$m, lv.formula = lv, t.formula = F$t, data = d, id = id, q = 20), error = function(e) NULL),
                           warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(x)) return(list(ok = FALSE))
  s <- summary(x)$mean.table; k <- which(rownames(s) == "time:trt")
  if (!length(k)) return(list(ok = TRUE, has_term = FALSE))   # term dropped by mm() (collinear design), as in the study rows
  list(ok = TRUE, has_term = TRUE, est = unname(s[k, 1]), se = unname(s[k, 2]), p = unname(s[k, 4]), code = unname(x$control["convergence_code"]),
       niter = unname(x$control["n_iter"]), alpha = paste(sprintf("%.10g", x$alpha), collapse = ";"), alpha_names = paste(names(x$alpha), collapse = ";"),
       status = if (length(warns)) paste(warns, collapse = " | ") else NA_character_)
}
slope_min <- function(a, nm) { v <- as.numeric(strsplit(a, ";")[[1]]); n <- strsplit(nm, ";")[[1]]
  k <- grepl("^log\\(sigma\\):", n) & !grepl("\\(Intercept\\)", n); if (!any(k)) NA_real_ else min(exp(v[k])) }
run_task <- function(r) {
  tk <- TASKS[r]; F <- FORM[[as.character(tk$model)]]; acc <- vector("list", CHUNK)
  for (j in seq_len(CHUNK)) {
    i <- (tk$chunk - 1L) * CHUNK + j
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 1000L * tk$model + i)
    d <- tryCatch(as.data.frame(make_all_data(tk$N, tk$ntime, i, tk$model, tk$label)), error = function(e) NULL)
    if (is.null(d)) next
    A <- fitw(d, F, F$lv)
    if (!A$ok) { acc[[j]] <- data.table(sim_id = i, A_ok = FALSE); next }
    if (!A$has_term) { acc[[j]] <- data.table(sim_id = i, A_ok = TRUE, A_has_term = FALSE); next }
    sm <- slope_min(A$alpha, A$alpha_names)
    bnd <- !(A$code %in% c(1L, 2L)) || (!is.na(sm) && sm < 1e-3)
    B <- if (bnd) fitw(d, F, ~1) else list(ok = NA)
    okB <- isTRUE(B$ok) && isTRUE(B$has_term)
    useB <- bnd && okB && is.finite(B$est) && is.finite(B$se)
    acc[[j]] <- data.table(sim_id = i, A_ok = TRUE, A_has_term = TRUE, A_est = A$est, A_se = A$se, A_p = A$p, A_code = A$code, A_niter = A$niter, A_alpha = A$alpha,
      A_alpha_names = A$alpha_names, A_status = A$status, A_slope_min = sm, bnd = bnd, B_ok = B$ok,
      B_est = if (okB) B$est else NA_real_, B_se = if (okB) B$se else NA_real_, B_p = if (okB) B$p else NA_real_,
      B_code = if (okB) B$code else NA_integer_, B_status = if (okB) B$status else NA_character_, used = if (useB) "B" else "A",
      cc_est = if (useB) B$est else A$est, cc_se = if (useB) B$se else A$se, cc_p = if (useB) B$p else A$p,
      cc_status = if (useB) B$status else A$status)
  }
  res <- rbindlist(Filter(Negate(is.null), acc), fill = TRUE)
  res[, `:=`(dgp = DGP, model = tk$model, N = tk$N, ntime = tk$ntime, label = tk$label, chunk = tk$chunk)]
  tmp <- paste0(tk$file, ".tmp"); saveRDS(res, tmp); file.rename(tmp, tk$file)
  nrow(res)
}
t0 <- Sys.time()
n <- parallel::mclapply(seq_len(nrow(TASKS)), function(r) tryCatch(run_task(r), error = function(e) { message("TASK FAILED ", basename(TASKS$file[r]), ": ", conditionMessage(e)); -1L }), mc.cores = CORES, mc.preschedule = FALSE)
cat(DGP, "done:", sum(unlist(n) >= 0), "tasks ok,", sum(unlist(n) < 0), "failed;", format(difftime(Sys.time(), t0, units = "hours")), "\n")
