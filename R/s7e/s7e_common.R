# =============================================================================
# s7e_common.R -- S7e: fourth data-generating family (Gaussian copula with logistic margins and serial latent
# correlation, R/v4/rundir/copula/) on the main grid (M1-M5 x 16 N x 6 T = 480 design points), H0 and H1,
# 1,000 replicates per arm, chunks of 100. One replicate:
#   1. degenerate-design flag (D-11; replays the design draws, identical to the other families);
#   2. seed of the study, 20260805 + 1000 * model + i, L'Ecuyer-CMRG (D-05, D-07): same designs as the other families;
#   3. the study's analyses through R/v4/rundir/do_analysis.R: mmm, gee (exchangeable, geepack), glmer (Wald, + shrink),
#      glmer_lrt (bobyqa, random structure of the Wald, D-08). The geesmv corrections are not run through gee_ss:
#      Mancl-DeRouen, Kauermann-Carroll and Fay-Graubard come from fewclust on the exchangeable fit (D-51; S7d showed
#      fewclust = geesmv for MD and FG), geesmv only on replicates 1-50 as cross-check (as D-48);
#   4. MMM with the convergence check of D-32 in M3-M5 (fit A = study fit, fit B with lv.formula = ~1 at the boundary);
#   5. the S7d strategies: geeglm with exchangeable, independence and AR(1) working correlation, fewclust md/kc/fg,
#      CR2 (exchangeable, independence), geesmv on replicates 1-50 (R/s7d/s7d_common.R::one_cor);
#   6. the sign-flip score test after set.seed(20261022 + i), primary and strict versions (D-46, D-47).
# =============================================================================
source(path.expand("~/MMM/R/s7d/s7d_common.R"))   # s3_common (libraries, seeds, design_degenerate), one_cor, s7d_flip, fewclust
COP_DIR <- file.path(RUNDIR_V4, "copula")
DGP_COP <- "Copula-AR1"
TABDIR <- file.path(REPO, "tables")
OUT_S7E <- file.path(REPO, "results_v4/S7e"); dir.create(file.path(OUT_S7E, "rows"), recursive = TRUE, showWarnings = FALSE)
TECHS_S7E <- c("mmm", "gee", "glmer", "glmer_lrt")
FORM_CC <- list(`3` = list(m = binY ~ time * trt, lv = ~time, t = ~1),
                `4` = list(m = binY ~ time * trt, lv = ~time, t = ~trt),
                `5` = list(m = binY ~ time * trt + dose, lv = ~time + dose, t = ~trt + steroids))

load_env_cop <- function() {
  env <- new.env(parent = globalenv())
  for (f in c("copula_params.R", "gen_copula_data.R", "generation_models_copula.R")) sys.source(file.path(COP_DIR, f), envir = env)
  for (f in c("manual/generate_data.R", "ar1/make_all_data.R")) sys.source(file.path(RUNDIR_V4, f), envir = env)
  for (f in ANALYSIS_FILES) sys.source(file.path(RUNDIR_V4, f), envir = env)
  env
}
# MMM fit with captured warnings (as R/s6/s6_mmm_cc_run.R::fitw)
fitw_cc <- function(d, F, lv) {
  warns <- character(0)
  x <- withCallingHandlers(tryCatch(mm(F$m, lv.formula = lv, t.formula = F$t, data = d, id = id, q = 20), error = function(e) NULL),
                           warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(x)) return(list(ok = FALSE))
  s <- summary(x)$mean.table; k <- which(rownames(s) == "time:trt")
  if (!length(k)) return(list(ok = TRUE, has_term = FALSE))
  list(ok = TRUE, has_term = TRUE, est = unname(s[k, 1]), se = unname(s[k, 2]), p = unname(s[k, 4]), code = unname(x$control["convergence_code"]),
       alpha = paste(sprintf("%.10g", x$alpha), collapse = ";"), alpha_names = paste(names(x$alpha), collapse = ";"),
       status = if (length(warns)) paste(warns, collapse = " | ") else NA_character_)
}
cc_fits_s7e <- function(d, model, i) {   # D-32 rule; called only for models 3-5
  F <- FORM_CC[[as.character(model)]]; A <- fitw_cc(d, F, F$lv)
  if (!A$ok) return(data.table(sim_id = i, A_ok = FALSE))
  if (!A$has_term) return(data.table(sim_id = i, A_ok = TRUE, A_has_term = FALSE))
  sm <- slope_min(A$alpha, A$alpha_names); bnd <- !(A$code %in% c(1L, 2L)) || (!is.na(sm) && sm < 1e-3)
  B <- if (bnd) fitw_cc(d, F, ~1) else list(ok = NA); okB <- isTRUE(B$ok) && isTRUE(B$has_term)
  useB <- bnd && okB && is.finite(B$est) && is.finite(B$se)
  data.table(sim_id = i, A_ok = TRUE, A_has_term = TRUE, A_est = A$est, A_se = A$se, A_p = A$p, A_code = A$code, A_status = A$status,
             A_slope_min = sm, bnd = bnd, B_ok = B$ok, B_est = if (okB) B$est else NA_real_, B_se = if (okB) B$se else NA_real_,
             B_p = if (okB) B$p else NA_real_, B_code = if (okB) B$code else NA_integer_, used = if (useB) "B" else "A",
             cc_est = if (useB) B$est else A$est, cc_se = if (useB) B$se else A$se, cc_p = if (useB) B$p else A$p,
             cc_status = if (useB) B$status else A$status)
}
s7e_tasks <- function() {
  x <- CJ(scen = sprintf("M%d", 1:5), N = as.integer(N_LEVELS), ntime = as.integer(T_LEVELS), label = c("mmlb", "manual"), chunk = 1:10)
  stopifnot(nrow(x) == 9600L)
  x[, model := as.integer(sub("M", "", scen))]
  x[, file := file.path(OUT_S7E, "rows", sprintf("copula_%s_N%d_T%d_%s_c%02d.rds", tolower(scen), N, ntime, label, chunk))]
  x
}
one_rep_s7e <- function(env, t, i, timing = FALSE) {
  tm <- c(gen = NA, study = NA, cc = NA, gee3 = NA, flip = NA)
  deg <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  t0 <- proc.time()[[3]]
  d <- tryCatch(as.data.frame(env$make_all_data(t$N, t$ntime, i, t$model, t$label)), error = function(e) e)
  tm["gen"] <- proc.time()[[3]] - t0
  rec <- data.table(sim_id = i, generated = !inherits(d, "error"), design_degenerate = deg,
                    gen_error = if (inherits(d, "error")) substr(conditionMessage(d), 1, 180) else NA_character_,
                    n_obs = NA_integer_, sum_y = NA_real_, p_flip = NA_real_, p_flip_strict = NA_real_, flip_why = NA_character_,
                    flip_strict_why = NA_character_, flip_max_abs_coef = NA_real_, flip_n_ref = NA_integer_)
  if (inherits(d, "error")) return(list(rec = rec, study = NULL, cc = NULL, gee = NULL, tm = tm))
  d <- d[order(d$id, d$time), ]; rec$n_obs <- nrow(d); rec$sum_y <- sum(d$binY)
  t0 <- proc.time()[[3]]
  st <- lapply(TECHS_S7E, function(tc) tryCatch(env$do_analysis(d, t$N, t$ntime, t$label, t$model, tc, i),
          error = function(e) trace_row(tc, t$N, t$ntime, t$label, t$model, i, paste("do_analysis error:", conditionMessage(e)))))
  st <- as.data.table(dplyr::bind_rows(st)); st <- st[is.na(term) | term == "time:trt"]
  tm["study"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
  cc <- if (t$model >= 3L) cc_fits_s7e(d, t$model, i) else NULL
  tm["cc"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
  form <- FORM_OF(t$scen)
  gee <- rbindlist(lapply(names(CORS), function(cn) one_cor(env, d, form, cn, i))); gee[, sim_id := i]
  tm["gee3"] <- proc.time()[[3]] - t0; t0 <- proc.time()[[3]]
  set.seed(FLIP_SEED(i))
  fl <- tryCatch(s7d_flip(d, t$scen), error = function(e) NULL)
  if (is.null(fl)) rec$flip_why <- "flip error" else {
    rec$p_flip <- fl$p[1]; rec$p_flip_strict <- fl$p[2]; rec$flip_why <- fl$why[1]; rec$flip_strict_why <- fl$why[2]
    rec$flip_max_abs_coef <- fl$max_abs_coef[1]; rec$flip_n_ref <- fl$n_ref[1] }
  tm["flip"] <- proc.time()[[3]] - t0
  list(rec = rec, study = st, cc = cc, gee = gee, tm = tm)
}
run_chunk_s7e <- function(t, ids = (t$chunk - 1L) * 100L + seq_len(100L), out_file = t$file) {
  env <- load_env_cop()
  out <- lapply(ids, function(i) one_rep_s7e(env, t, i))
  res <- lapply(c("rec", "study", "cc", "gee"), function(k) rbindlist(lapply(out, `[[`, k), fill = TRUE)); names(res) <- c("recs", "study", "cc", "gee")
  for (k in names(res)) if (nrow(res[[k]])) res[[k]][, `:=`(dgp = DGP_COP, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
  res$secs <- do.call(rbind, lapply(out, `[[`, "tm"))
  if (!is.null(out_file)) { tmp <- paste0(out_file, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out_file) }
  invisible(res)
}
