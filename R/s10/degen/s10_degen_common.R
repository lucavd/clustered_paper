# =============================================================================
# s10_degen_common.R -- S10, track "Degeneri": regeneration of every dataset of the S7d run (1,896 points) and of the
# S7e run (480 Copula-AR1 points), 2 arms x 1,000 replicates, exactly as one_rep_s7d / one_rep_s7e
# (design_degenerate() first, then set_rep_seed(model, i), then the generator), and per dataset:
#   * GD2 inputs: generated, design_degenerate, n_obs, sum_y, next to the same fields of the stored recs;
#   * pre-registered classes (docs/steps/S10.md, "Attese pre-registrate"), on generated, non-degenerate datasets
#     without constant outcome: S (LP of Konis 2007, function lp_sep copied verbatim from R/s10/probe/s10_probe2.R,
#     formula FORM_OF(scen), independence), Si (detectseparation 0.4.0 declares time:trt infinite; run on every
#     LP-separated dataset and wherever LP returns no solution), V1, E1, B0, B1; class L is per strategy (s10_degen_join.R);
#   * descriptives: occasions without events / with only events; arm x occasion cells empty, without events, only events.
# Read-only on everything outside results_v4/S10/degen. One output file per chunk of 100 replicates.
# =============================================================================
suppressMessages(source(path.expand("~/MMM/R/s7e/s7e_common.R")))   # sources s7d_common.R (s3_common, FORM_OF, s7d_tasks)
.libPaths(c(.libPaths(), path.expand("~/R/s10_lib")))               # after s7d_common.R, which resets the paths
suppressMessages({ library(detectseparation); library(lpSolveAPI) })
stopifnot(as.character(packageVersion("detectseparation")) == "0.4.0")
data.table::setDTthreads(1L)
OUT_DG <- file.path(REPO, "results_v4/S10/degen")
for (p in c("rows", "logs", "gd3")) dir.create(file.path(OUT_DG, p), recursive = TRUE, showWarnings = FALSE)
# ---- lp_sep: verbatim from the probe (extracted from the file, so it cannot drift) -------------------------------
local({ src <- readLines(file.path(REPO, "R/s10/probe/s10_probe2.R")); a <- grep("^lp_sep <- function", src)
  b <- a + which(src[(a + 1):length(src)] == "}")[1]; eval(parse(text = src[a:b]), envir = globalenv()) })
stopifnot(exists("lp_sep"))
dg_tasks <- function() {
  a <- s7d_tasks()[, .(src = "S7d", dgp, scen, N, ntime, label, chunk, model, rec_file = file)]
  b <- s7e_tasks()[, .(src = "S7e", dgp = DGP_COP, scen, N, ntime, label, chunk, model, rec_file = file)]
  x <- rbind(a, b); x[, out := file.path(OUT_DG, "rows", basename(rec_file))]
  stopifnot(nrow(x) == 47520L, !anyDuplicated(x$out)); x
}
ENVS <- new.env()
dg_env <- function(dgp) { if (is.null(ENVS[[dgp]])) ENVS[[dgp]] <- if (dgp == DGP_COP) load_env_cop() else load_env(dgp, RUNDIR_V4); ENVS[[dgp]] }
dg_gen <- function(t) { env <- dg_env(t$dgp)
  if (t$src == "S7e") function(N, ntime, i, label) env$make_all_data(N, ntime, i, t$model, label) else s7d_gen(env, t) }
BLANK <- list(sim_id = NA_integer_, generated = NA, design_degenerate = NA, gen_error = NA_character_, n_obs = NA_integer_, sum_y = NA_real_,
  const = NA, elig = NA, time_ok = NA, n_subj = NA_integer_, n_ev_subj = NA_integer_, n_ne_subj = NA_integer_, n_var_subj = NA_integer_,
  ev_arm0 = NA_integer_, ev_arm1 = NA_integer_, nobs_arm0 = NA_integer_, nobs_arm1 = NA_integer_, m2_arm0 = NA_integer_, m2_arm1 = NA_integer_,
  occ_noev = NA_integer_, occ_allev = NA_integer_, cell_empty = NA_integer_, cell_noev = NA_integer_, cell_allev = NA_integer_,
  V1 = NA, E1 = NA, B0 = NA, B1 = NA, lp = NA, S = NA, ds_run = FALSE, ds_outcome = NA, ds_err = NA_character_, Si = NA,
  ds_inf_terms = NA_character_, t_gen = NA_real_, t_lp = NA_real_, t_ds = NA_real_)
# detectseparation on one dataset: outcome and infinite coefficients
dg_ds <- function(form, d) {
  s <- tryCatch(suppressWarnings(glm(form, data = d, family = binomial, method = "detect_separation")), error = function(e) e)
  if (inherits(s, "error")) return(list(outcome = NA, err = substr(conditionMessage(s), 1, 160), Si = NA, inf = NA_character_))
  cf <- s$coefficients; k <- match("time:trt", names(cf))
  list(outcome = isTRUE(s$outcome), err = NA_character_, Si = if (is.na(k)) NA else is.infinite(cf[[k]]),
       inf = paste(names(cf)[is.infinite(cf)], collapse = ";"))
}
dg_one <- function(gen, t, i, ds_always = FALSE) {
  r <- BLANK; r$sim_id <- as.integer(i)
  r$design_degenerate <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  t0 <- proc.time()[[3]]
  d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) e)
  r$t_gen <- proc.time()[[3]] - t0
  r$generated <- !inherits(d, "error")
  if (!r$generated) { r$gen_error <- substr(conditionMessage(d), 1, 180); return(r) }
  d <- d[order(d$id, d$time), ]; y <- d$binY
  r$n_obs <- nrow(d); r$sum_y <- sum(y)
  r$const <- r$sum_y == 0 || r$sum_y == r$n_obs
  r$elig <- !r$design_degenerate && !r$const
  # subject-level descriptives
  id <- factor(d$id); ev <- tapply(y, id, max); ne <- tapply(1 - y, id, max); ns <- tapply(y, id, length)
  arm_s <- tapply(d$trt, id, function(z) z[1])
  r$n_subj <- length(ev); r$n_ev_subj <- as.integer(sum(ev)); r$n_ne_subj <- as.integer(sum(ne)); r$n_var_subj <- as.integer(sum(ev & ne))
  af <- factor(d$trt, levels = 0:1)
  ea <- tapply(y, af, sum); na <- tapply(y, af, length); ea[is.na(ea)] <- 0L; na[is.na(na)] <- 0L
  m2 <- tapply(ns >= 2, factor(arm_s, levels = 0:1), sum); m2[is.na(m2)] <- 0L
  r$ev_arm0 <- as.integer(ea[["0"]]); r$ev_arm1 <- as.integer(ea[["1"]]); r$nobs_arm0 <- as.integer(na[["0"]]); r$nobs_arm1 <- as.integer(na[["1"]])
  r$m2_arm0 <- as.integer(m2[["0"]]); r$m2_arm1 <- as.integer(m2[["1"]])
  occ <- seq.int(0L, t$ntime - 1L); r$time_ok <- all(d$time %in% occ)
  tf <- factor(d$time, levels = occ)
  eo <- tapply(y, tf, sum); no <- tapply(y, tf, length); eo[is.na(eo)] <- 0; no[is.na(no)] <- 0
  r$occ_noev <- as.integer(sum(no > 0 & eo == 0)); r$occ_allev <- as.integer(sum(no > 0 & eo == no))
  ec <- tapply(y, list(af, tf), sum); nc <- tapply(y, list(af, tf), length); ec[is.na(ec)] <- 0; nc[is.na(nc)] <- 0
  r$cell_empty <- as.integer(sum(nc == 0)); r$cell_noev <- as.integer(sum(nc > 0 & ec == 0)); r$cell_allev <- as.integer(sum(nc > 0 & ec == nc))
  # pre-registered classes (defined on eligible datasets only)
  if (r$elig) {
    r$V1 <- r$n_var_subj <= 1L; r$E1 <- min(r$n_ev_subj, r$n_ne_subj) <= 1L
    r$B0 <- any(ea == 0 | ea == na); r$B1 <- any(m2 == 1L)
  }
  if (r$elig || ds_always) {
    form <- FORM_OF(t$scen); X <- model.matrix(form, d)
    t0 <- proc.time()[[3]]; a <- lp_sep(X, y); r$t_lp <- proc.time()[[3]] - t0
    r$lp <- a; r$S <- if (r$elig) isTRUE(a) else NA
    if (isTRUE(a) || is.na(a) || ds_always) {
      t0 <- proc.time()[[3]]; s <- dg_ds(form, d); r$t_ds <- proc.time()[[3]] - t0
      r$ds_run <- TRUE; r$ds_outcome <- s$outcome; r$ds_err <- s$err; r$ds_inf_terms <- s$inf
      if (r$elig) r$Si <- isTRUE(a) && isTRUE(s$Si)
    } else if (r$elig) r$Si <- FALSE
  }
  r
}
dg_chunk <- function(t) {
  t0 <- proc.time()[[3]]; gen <- dg_gen(t)
  ids <- (t$chunk - 1L) * 100L + seq_len(100L)
  X <- rbindlist(lapply(ids, function(i) dg_one(gen, t, i)))
  R <- readRDS(t$rec_file)$recs[, .(sim_id = as.integer(sim_id), rec_generated = generated, rec_deg = design_degenerate,
                                     rec_n_obs = as.integer(n_obs), rec_sum_y = as.numeric(sum_y))]
  X <- merge(X, R, by = "sim_id", all = TRUE, sort = TRUE)
  X[, `:=`(src = t$src, dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
  attr(X, "secs") <- proc.time()[[3]] - t0
  tmp <- paste0(t$out, ".tmp"); saveRDS(X, tmp); file.rename(tmp, t$out)
  invisible(nrow(X))
}
