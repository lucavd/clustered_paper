# =============================================================================
# s7g_common.R -- S7g: the study's mixed model (glmer, random intercept only: M1, M2, M1c) refitted with adaptive
# Gauss-Hermite quadrature on the study's own datasets, H0 only (label "manual").
#   Design points: M1, M2 x 4 mechanisms (GenBinaryY, GLMM-native, GLMM-AR1, Copula-AR1) x N {10, 25, 55, 145} x T {2, 4, 10};
#                  M1c x 3 mechanisms (M1c does not exist for Copula-AR1) x N {10, 25, 45} x T {2, 4, 6}
#                  -> 123 design points, 1,000 replicates each, chunks of 100 (1,230 chunk files).
#   Per replicate: degenerate-design flag (D-11); study seed 20260805 + 1000 * model + i with L'Ecuyer-CMRG (D-05, D-07;
#                  M1c uses the seeds of M1); dataset regenerated as in S3 (three study mechanisms) and as in S7e
#                  (Copula-AR1, sorted by id and time as there); then the study's glmer through R/v4/rundir/do_analysis.R
#                  with nAGQ = 1 (Laplace: must reproduce the stored rows), nAGQ = 10 and, for replicates 1-100, nAGQ = 25.
#   Only nAGQ changes: nAGQ is injected by masking `glmer` in the family environment (the enclosure of the analysis
#   functions), so formula, data, family and glmerControl are those of the study.
# =============================================================================
source(path.expand("~/MMM/R/s7e/s7e_common.R"))   # -> s7d_common -> s3_common: load_env, gen_data, design_degenerate, set_rep_seed
OUT_S7G <- file.path(REPO, "results_v4/S7g"); dir.create(file.path(OUT_S7G, "rows"), recursive = TRUE, showWarnings = FALSE)
N_G <- c(10L, 25L, 55L, 145L); T_G <- c(2L, 4L, 10L); N_G1C <- c(10L, 25L, 45L); T_G1C <- c(2L, 4L, 6L)
NAGQ_MAIN <- 10L; NAGQ_CHECK <- 25L; CHECK_MAX_I <- 100L
s7g_points <- function() {
  a <- CJ(dgp = c(DGPS, DGP_COP), scen = c("M1", "M2"), N = N_G, ntime = T_G)
  b <- CJ(dgp = DGPS, scen = "M1c", N = N_G1C, ntime = T_G1C)
  x <- rbind(a, b); stopifnot(nrow(x) == 123L); x
}
s7g_tasks <- function() {
  x <- s7g_points()[, .(chunk = 1:10), by = .(dgp, scen, N, ntime)]
  x[, `:=`(label = "manual", model = as.integer(sub("M", "", sub("M1c", "M1", scen, fixed = TRUE))))]
  x[, file := file.path(OUT_S7G, "rows", sprintf("%s_%s_N%d_T%d_manual_c%02d.rds", dgp_tag(dgp), tolower(scen), N, ntime, chunk))]
  stopifnot(nrow(x) == 1230L); x
}
s7g_env <- function(dgp) if (dgp == DGP_COP) load_env_cop() else load_env(dgp, RUNDIR_V4)
fit_glmer_k <- function(env, d, t, i, K) {
  stopifnot(!exists("glmer", envir = env, inherits = FALSE))
  if (K > 1L) { assign("glmer", function(...) lme4::glmer(..., nAGQ = K), envir = env); on.exit(rm("glmer", envir = env)) }
  t0 <- proc.time()[[3]]
  r <- tryCatch(env$do_analysis(d, t$N, t$ntime, t$label, t$model, "glmer", i),
                error = function(e) trace_row("glmer", t$N, t$ntime, t$label, t$model, i, paste("do_analysis error:", conditionMessage(e))))
  secs <- proc.time()[[3]] - t0
  r <- as.data.table(r)[method == "glmer" & (is.na(term) | term == "time:trt")]
  if (!nrow(r)) r <- data.table(term = NA_character_, estimate = NA_real_, se = NA_real_, p = NA_real_, converged = FALSE, status = "no time:trt row")
  if (!"D" %in% names(r)) r[, D := NA_real_]
  r[, .(sim_id = i, nAGQ = K, term, estimate, se, p, D, converged, status, secs = secs)]
}
one_rep_s7g <- function(env, t, i) {
  deg <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  d <- tryCatch(as.data.frame(gen_data(env, t$scen, t$model, t$N, t$ntime, i, t$label)), error = function(e) e)
  rec <- data.table(sim_id = i, generated = !inherits(d, "error"), design_degenerate = deg, n_obs = NA_integer_, sum_y = NA_real_)
  if (inherits(d, "error")) return(list(rec = rec, rows = NULL))
  if (t$dgp == DGP_COP) d <- d[order(d$id, d$time), ]
  rec$n_obs <- nrow(d); rec$sum_y <- sum(d$binY)
  Ks <- c(1L, NAGQ_MAIN, if (i <= CHECK_MAX_I) NAGQ_CHECK)
  list(rec = rec, rows = rbindlist(lapply(Ks, function(K) fit_glmer_k(env, d, t, i, K)), fill = TRUE))
}
run_chunk_s7g <- function(t, ids = (t$chunk - 1L) * 100L + seq_len(100L), out_file = t$file) {
  env <- s7g_env(t$dgp)
  out <- lapply(ids, function(i) one_rep_s7g(env, t, i))
  res <- list(recs = rbindlist(lapply(out, `[[`, "rec"), fill = TRUE), rows = rbindlist(lapply(out, `[[`, "rows"), fill = TRUE))
  for (k in names(res)) if (nrow(res[[k]])) res[[k]][, `:=`(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
  if (!is.null(out_file)) { tmp <- paste0(out_file, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out_file) }
  invisible(res)
}
