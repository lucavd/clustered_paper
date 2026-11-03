# =============================================================================
# s7c_common.R -- shared setup for S7c (pilot of tests whose null is the paper's).
# Datasets are regenerated exactly as R/s3/s3_common.R::one_rep (degenerate-design flag first, then seeded generation).
# =============================================================================
source(path.expand("~/MMM/R/s3/s3_common.R"))   # libraries, load_env, set_rep_seed, design_degenerate, gen_data, dgp_tag
source(path.expand("~/MMM/R/s7c/s7c_tests.R"))
suppressMessages(library(data.table))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
data.table::setDTthreads(1L)
OUT_S7C <- file.path(REPO, "results_v4/S7c"); dir.create(file.path(OUT_S7C, "pilot"), recursive = TRUE, showWarnings = FALSE)
S7C_LIB <- path.expand("~/R/s7c_lib")            # flipscores, fewclust (kept out of the project library)
RESAMPLE_SEED <- function(i) 20261022L + as.integer(i)
PILOT_N <- c(10L, 25L); PILOT_T <- c(2L, 3L, 10L)
pilot_tasks <- function() {
  x <- CJ(dgp = DGPS, N = PILOT_N, ntime = PILOT_T, label = c("manual", "mmlb"), chunk = 1:10)
  x[, file := file.path(OUT_S7C, "pilot", sprintf("%s_m1_N%d_T%d_%s_c%02d.rds", dgp_tag(dgp), N, ntime, label, chunk))]
  x[, orig := file.path(REPO, "results_full", sprintf("%s_m1_N%d_T%d_%s_c%02d.rds", dgp_tag(dgp), N, ntime, label, chunk))]
  x
}
gee_tt <- function(env, d, N, ntime, label) {
  g <- tryCatch(env$analyze_model1_gee(d, N, ntime, label), error = function(e) NULL)
  if (is.null(g)) return(c(est = NA_real_, se = NA_real_))
  k <- which(g$term == "time:trt"); if (!length(k)) return(c(est = NA_real_, se = NA_real_))
  c(est = g$estimate[k], se = g$se[k])
}
one_rep_s7c <- function(env, N, ntime, label, i) {
  deg <- design_degenerate(1L, N, ntime, i)
  set_rep_seed(1L, i)
  d <- tryCatch(as.data.frame(gen_data(env, "M1", 1L, N, ntime, i, label)), error = function(e) e)
  if (inherits(d, "error")) return(list(rec = data.table(sim_id = i, generated = FALSE, design_degenerate = deg,
                                                         gen_error = substr(conditionMessage(d), 1, 180)), tests = NULL))
  gt <- gee_tt(env, d, N, ntime, label)
  set.seed(RESAMPLE_SEED(i))
  tt <- tryCatch(run_tests(d), error = function(e) data.frame(test = c("flip", "subj_t", "subj_perm", "naive_perm", "naive_slope"),
                 p = NA_real_, why = paste("run_tests error:", substr(conditionMessage(e), 1, 150))))
  tt <- as.data.table(tt); tt[, sim_id := i]
  list(rec = data.table(sim_id = i, generated = TRUE, design_degenerate = deg, gen_error = NA_character_,
                        n_obs = nrow(d), sum_y = sum(d$binY), gee_est = gt[["est"]], gee_se = gt[["se"]]),
       tests = tt)
}
run_chunk_s7c <- function(t) {
  env <- load_env(t$dgp, RUNDIR_V4)
  ids <- (t$chunk - 1L) * 100L + seq_len(100L)
  out <- lapply(ids, function(i) one_rep_s7c(env, t$N, t$ntime, t$label, i))
  recs <- rbindlist(lapply(out, `[[`, "rec"), fill = TRUE); tests <- rbindlist(lapply(out, `[[`, "tests"), fill = TRUE)
  for (D in list(recs, tests)) D[, `:=`(dgp = t$dgp, N = t$N, ntime = t$ntime, label = t$label)]
  tmp <- paste0(t$file, ".tmp"); saveRDS(list(recs = recs, tests = tests), tmp); file.rename(tmp, t$file)
  nrow(tests)
}
