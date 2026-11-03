# =============================================================================
# s3_common.R -- functions shared by the S3 driver and its checks.
# One replicate is produced exactly as run_full.R::run_chunk() does
# (set.seed -> generate -> analyses, nothing in between), with three v4 additions:
# explicit RNGkind (D-07), a trace row for every failure (E3), and per-replicate
# records (design_degenerate for D-11, design identity against a reference
# generator for Q4). Reference generation and the design flag are computed
# BEFORE the seeded generation of the analysed dataset, so the random-number
# sequence seen by the analyses is the one of run_full.R.
# =============================================================================
REPO   <- path.expand("~/MMM")
LIB    <- file.path(REPO, "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({
  library(MMLB); library(lme4); library(geepack); library(geesmv)
  library(dplyr); library(tibble); library(MASS); library(fastGHQuad)
})
source(file.path(REPO, "R/protocol/dgm_common.R"))     # DGPS, FAMILY_FILES, seed_of, load_family
RUNDIR_ORIG <- file.path(REPO, ".rundir")
RUNDIR_V4   <- file.path(REPO, "R/v4/rundir")
OUT_S3      <- file.path(REPO, "results_v4/S3")
RNG_PROD    <- c("L'Ecuyer-CMRG", "Inversion", "Rejection")   # D-07
TECHS_ALL   <- c("mmm", "gee", "glmer", "gee_ss", "glmer_lrt")
TECHS_LRT   <- c("glmer_lrt", "gee")                          # D-15: gee = data-identity control
ANALYSIS_FILES <- c(sprintf("analysis_model%d.R", 1:5), "analysis_smallsample.R", "do_analysis.R")
M1C_FILE <- c("GenBinaryY" = "generation_model1c.R", "GLMM-native" = "manual/generation_model1c.R",
              "GLMM-AR1" = "ar1/generation_model1c.R")
N_LEVELS <- c(10, 15, 20, seq(25, 145, 10)); T_LEVELS <- c(2, 3, 4, 6, 8, 10)
N_M1C <- c(10, 15, 20, 25, 35, 45); T_M1C <- c(2, 3, 4, 6)

set_rep_seed <- function(model, i) { RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(seed_of(model, i)) }

load_env <- function(dgp, rundir) {
  env <- load_family(dgp, rundir)
  if (file.exists(file.path(rundir, "m1c_params.R"))) {
    sys.source(file.path(rundir, "m1c_params.R"), envir = env)
    sys.source(file.path(rundir, M1C_FILE[[dgp]]), envir = env)
  }
  for (f in ANALYSIS_FILES) sys.source(file.path(rundir, f), envir = env)
  env
}

gen_data <- function(env, scen, model, N, ntime, i, label) {
  if (scen == "M1c") env$generate_data_model1c(N, ntime, i, label)
  else env$make_all_data(N, ntime, i, model, label)
}

# D-11: estimable design = both arms observed beyond time 0. Replays the first two
# design draws of every generator (nclust, trt), which are identical across
# scenarios, mechanisms and arms (S1).
design_degenerate <- function(model, N, ntime, i) {
  set_rep_seed(model, i)
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))
  trt <- rbinom(N, 1, 0.5)
  !(any(trt == 1 & nclust > 1) && any(trt == 0 & nclust > 1))
}

trace_row <- function(method, N, ntime, label, model, i, msg) tibble::tibble(
  term = NA_character_, estimate = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_,
  method = method, N = N, ntime = ntime, label = label, model = model, sim_id = i,
  converged = FALSE, status = paste("[v4 trace]", substr(msg, 1, 180)))

# ref: NULL, or list(env, scen, model, label) of a reference generator whose
# design columns must coincide (Q4).
one_rep <- function(env, scen, model, N, ntime, label, i, techs, ref = NULL) {
  deg <- design_degenerate(model, N, ntime, i)
  dref <- NULL
  if (!is.null(ref)) {
    set_rep_seed(ref$model, i)
    dref <- tryCatch(gen_data(ref$env, ref$scen, ref$model, N, ntime, i, ref$label), error = function(e) NULL)
  }
  set_rep_seed(model, i)
  d <- tryCatch(gen_data(env, scen, model, N, ntime, i, label), error = function(e) e)
  rec <- data.frame(sim_id = i, generated = !inherits(d, "error"),
                    gen_error = if (inherits(d, "error")) substr(conditionMessage(d), 1, 180) else NA_character_,
                    design_degenerate = deg, design_identical_ref = NA,
                    calib_max_error = NA_real_, n_obs = NA_integer_, sum_y = NA_real_)
  if (inherits(d, "error"))
    return(list(rows = trace_row("generation", N, ntime, label, model, i, paste("generation error:", conditionMessage(d))), rec = rec))
  rec$n_obs <- nrow(d); rec$sum_y <- sum(d$binY)
  if (!is.null(attr(d, "calib_max_error"))) rec$calib_max_error <- attr(d, "calib_max_error")
  if (!is.null(ref)) {
    cols <- intersect(c("id", "time", "trt", "dose"), names(d))
    rec$design_identical_ref <- !is.null(dref) &&
      identical(lapply(as.data.frame(d)[cols], as.vector), lapply(as.data.frame(dref)[cols], as.vector))
  }
  rows <- lapply(techs, function(tc) tryCatch(env$do_analysis(d, N, ntime, label, model, tc, i),
    error = function(e) trace_row(tc, N, ntime, label, model, i, paste("do_analysis error:", conditionMessage(e)))))
  list(rows = dplyr::bind_rows(rows), rec = rec)
}

# ---- task table -------------------------------------------------------------
dgp_tag <- function(dgp) gsub("[^A-Za-z0-9]", "", dgp)
s3_tasks <- function() {
  g <- function(block, dgp, scen, model, N, Tt, labels, techs) {
    x <- expand.grid(block = block, dgp = dgp, scen = scen, N = N, ntime = Tt, label = labels,
                     chunk = 1:10, stringsAsFactors = FALSE)
    x$model <- ifelse(x$scen == "M1c", 1L, as.integer(sub("M", "", x$scen)))
    x$techs <- techs; x
  }
  M <- sprintf("M%d", 3:5)
  x <- rbind(
    g("F1", "GenBinaryY", M, NA, N_LEVELS, T_LEVELS, "manual", "all"),
    g("F2", "GLMM-AR1", "M5", NA, N_LEVELS, T_LEVELS, c("mmlb", "manual"), "all"),
    g("L", "GenBinaryY", M, NA, N_LEVELS, T_LEVELS, "mmlb", "lrt"),
    g("L", "GLMM-native", M, NA, N_LEVELS, T_LEVELS, c("mmlb", "manual"), "lrt"),
    g("L", "GLMM-AR1", c("M3", "M4"), NA, N_LEVELS, T_LEVELS, c("mmlb", "manual"), "lrt"),
    do.call(rbind, lapply(DGPS, function(dg) g("C", dg, "M1c", NA, N_M1C, T_M1C, c("mmlb", "manual"), "all"))))
  x$key  <- sprintf("%s_%s_N%d_T%d_%s_c%02d", dgp_tag(x$dgp), tolower(x$scen), x$N, x$ntime, x$label, x$chunk)
  x$file <- file.path(OUT_S3, x$block, paste0(x$key, ".rds"))
  x$orig_file <- ifelse(x$scen == "M1c", NA_character_,
    file.path(REPO, "results_full", sprintf("%s_m%d_N%d_T%d_%s_c%02d.rds", dgp_tag(x$dgp), x$model, x$N, x$ntime, x$label, x$chunk)))
  x
}

# ---- one chunk (runs in a worker) ---------------------------------------------
run_chunk_s3 <- function(t) {
  env  <- load_env(t$dgp, RUNDIR_V4)
  ref <- NULL
  if (t$block %in% c("F1", "F2")) ref <- list(env = load_env(t$dgp, RUNDIR_ORIG), scen = t$scen, model = t$model, label = t$label)
  if (t$block == "C") ref <- list(env = load_env(t$dgp, RUNDIR_ORIG), scen = "M1", model = 1L, label = t$label)
  techs <- if (t$techs == "all") TECHS_ALL else TECHS_LRT
  ids <- (t$chunk - 1L) * 100L + seq_len(100L)
  out <- lapply(ids, function(i) one_rep(env, t$scen, t$model, t$N, t$ntime, t$label, i, techs, ref))
  rows <- dplyr::bind_rows(lapply(out, `[[`, "rows"))
  recs <- do.call(rbind, lapply(out, `[[`, "rec"))
  rows$dgp <- t$dgp; rows$scenario <- t$scen; rows$block <- t$block; rows$gen_version <- "v4"
  rows$design_degenerate <- recs$design_degenerate[match(rows$sim_id, recs$sim_id)]
  rows$re_structure <- ifelse(rows$method == "glmer_lrt", env$lrt_re_structure(t$model), NA_character_)
  attr(rows, "sims") <- recs
  dir.create(dirname(t$file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(t$file, ".tmp"); saveRDS(rows, tmp, compress = "xz"); file.rename(tmp, t$file)
  c(rows = nrow(rows))
}
