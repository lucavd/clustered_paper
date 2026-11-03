# =============================================================================
# extract_dgm_params.R -- S1 of the Morris (2019) roadmap.
#
# Builds tables/V4_dgm_parameters.csv, the table of every data-generating
# parameter, READ FROM THE CODE THAT RAN rather than typed by hand:
#   * engine arguments are captured by wrapping GenBinaryY / gen_glmm_data /
#     gen_ar1_data / calibrate_beta_c_ar1 while make_all_data() runs;
#   * the design (follow-up truncation, covariate draws, their order) comes from
#     a static walk of the generator body, matched against the random-number
#     calls actually executed;
#   * the seed rule is parsed from run_full.R in the same staged directory.
# Capture is repeated at three cells; every specified parameter must be the
# same at all three (the parameters do not depend on N, T or the replicate).
#
# Usage: Rscript --vanilla R/protocol/extract_dgm_params.R [rundir] [out.csv]
# =============================================================================

LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(MMLB); library(dplyr); library(MASS); library(fastGHQuad)})

argv   <- commandArgs(trailingOnly = TRUE)
REPO   <- path.expand("~/MMM")
RUNDIR <- if (length(argv) >= 1) argv[1] else file.path(REPO, ".rundir")
OUT    <- if (length(argv) >= 2) argv[2] else file.path(REPO, "tables/V4_dgm_parameters.csv")
source(file.path(REPO, "R/protocol/dgm_common.R"))

CELLS <- data.frame(N = c(45, 10, 145), ntime = c(6, 3, 10), i = c(1, 500, 999))

# ---- seed rule, parsed from the staged driver --------------------------------
find_calls <- function(e, name, acc = list()) {
  if (is.call(e)) {
    if (identical(e[[1]], as.name(name))) acc[[length(acc) + 1L]] <- e
    for (k in seq_along(e)[-1]) if (!identical(e[[k]], quote(expr = ))) acc <- find_calls(e[[k]], name, acc)
  }
  acc
}
drv <- parse(file.path(RUNDIR, "run_full.R"))
ss  <- unlist(lapply(drv, find_calls, name = "set.seed"), recursive = FALSE)
stopifnot(length(ss) == 1L)
seed_expr <- ss[[1]][[2]]
seed_txt  <- paste(deparse(seed_expr), collapse = "")
for (m in 1:5) for (i in c(1, 1000))
  stopifnot(eval(seed_expr, list(model = m, i = i)) == seed_of(m, i))

# ---- RNG kind of the production workers -------------------------------------------
# run_full.R never calls RNGkind(); the replicates ran inside crew/mirai workers,
# which initialise L'Ecuyer-CMRG streams. The kind is therefore read from a real
# crew worker started exactly as run_full.R starts them, not from this session.
worker_rngkind <- function() {
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  suppressMessages(library(crew))
  ctrl <- crew::crew_controller_local(workers = 1, seconds_idle = 30)
  ctrl$start(); on.exit(ctrl$terminate())
  ctrl$push(command = paste(RNGkind(), collapse = ", "), name = "rngkind")
  ctrl$wait(mode = "all", seconds_timeout = 120)
  p <- ctrl$pop()
  stopifnot(!is.null(p), is.na(p$error[[1]]) || is.null(p$error[[1]]))
  p$result[[1]]
}
WORKER_RNGKIND <- worker_rngkind()
SESSION_RNGKIND <- paste(RNGkind(), collapse = ", ")
cat("RNGkind in a crew worker:", WORKER_RNGKIND, "| in this session:", SESSION_RNGKIND, "\n")
set_rngkind <- function(txt) {
  k <- trimws(strsplit(txt, ",")[[1]])
  RNGkind(kind = k[1], normal.kind = k[2], sample.kind = k[3])
}
set_rngkind(WORKER_RNGKIND)

# ---- one capture ----------------------------------------------------------------
capture_one <- function(dgp, model, arm, N, ntime, i) {
  env <- load_family(dgp, RUNDIR)
  rec <- install_recorders(env)
  set.seed(seed_of(model, i))
  out <- env$make_all_data(N, ntime, i, model, ARMS[[arm]])
  calls <- rec$engine_calls
  eng   <- Filter(function(cl) cl$engine == ENGINE_OF[[dgp]], calls)
  stopifnot(length(eng) == 1L)
  a <- eng[[1]]$args
  dat <- a$data

  # static design and its match with the executed draws
  fun <- if (dgp == "GLMM-AR1") env$.ar1_design else env[[paste0("generate_data_model", model)]]
  st  <- static_design(fun)
  dyn <- rec$rng
  pos <- 0L; matched <- integer(0)
  for (k in seq_along(dyn)) {
    cand <- which(vapply(st$draws, `[[`, "", "fn") == dyn[[k]]$fn)
    cand <- cand[cand > pos]
    stopifnot(length(cand) > 0L)
    pos <- cand[1]; matched <- c(matched, pos)
  }
  n_obs <- nrow(dat)
  rows <- list()
  add <- function(block, parameter, value, terms = "", note = "")
    rows[[length(rows) + 1L]] <<- data.frame(dgp = dgp, model = model, arm = arm,
      label = ARMS[[arm]], block = block, parameter = parameter, value = value,
      terms = terms, note = note, stringsAsFactors = FALSE)

  add("engine", "engine", ENGINE_OF[[dgp]])
  mf <- a$mean.formula
  X  <- model.matrix(mf, model.frame(mf, dat))
  add("mean", "mean.formula", f_txt(mf))
  if (dgp == "GenBinaryY") {
    add("mean", "beta", num_txt(a$beta), paste(colnames(X), collapse = ", "), "marginal")
    add("association", "lv.formula", f_txt(a$lv.formula))
    add("association", "t.formula",  f_txt(a$t.formula))
    lvX <- model.matrix(a$lv.formula, model.frame(a$lv.formula, dat))
    tX  <- model.matrix(a$t.formula,  model.frame(a$t.formula,  dat))
    add("association", "sigma", num_txt(a$sigma), paste(colnames(lvX), collapse = ", "),
        "linear predictor of the latent-variable loading: sigma_ij = x_lv' sigma (MMLB GenBinaryY)")
    add("association", "gamma", num_txt(a$gamma), paste(colnames(tX), collapse = ", "),
        "linear predictor of the transition (serial) coefficient: gamma_ij = x_t' gamma")
    add("association", "q", num_txt(a$q), "", "quadrature points in GenBinaryY")
  } else {
    D <- as.matrix(a$D)
    Z <- model.matrix(a$re.formula, model.frame(a$re.formula, dat))
    add("random_effects", "re.formula", f_txt(a$re.formula))
    add("random_effects", "D", mat_txt(D), paste(colnames(Z), collapse = ", "),
        "covariance of b_i, row-major")
    if (dgp == "GLMM-native") {
      add("mean", "beta_c", num_txt(a$beta_c), paste(colnames(X), collapse = ", "),
          "conditional, fixed constant")
    } else {
      cal <- Filter(function(cl) cl$engine == "calibrate_beta_c_ar1", calls)
      stopifnot(length(cal) == 1L)
      ca <- cal[[1]]$args
      stopifnot(isTRUE(all.equal(as.matrix(ca$D), D, tolerance = 0)),
                identical(ca$sigma2_ar, a$sigma2_ar))
      add("mean", "target_beta_m", num_txt(ca$target_beta_m),
          paste(colnames(X), collapse = ", "),
          "marginal target; beta_c is calibrated to it on every replicate")
      add("serial", "rho", num_txt(a$rho), "", "AR(1) correlation of latent residual")
      add("serial", "sigma2_ar", num_txt(a$sigma2_ar), "", "AR(1) latent residual variance")
      add("derived", "beta_c_calibrated", num_txt(a$beta_c), paste(colnames(X), collapse = ", "),
          sprintf("replicate-specific; this row at N=%d, T=%d, i=%d", N, ntime, i))
    }
  }

  # design
  add("design", "p_full_followup", num_txt(st$threshold), "nclust",
      "P(subject observed at all ntime occasions)")
  nc <- Filter(function(z) z$target == "nclust", st$assigns)
  add("design", "followup_rule", nc[[length(nc)]]$rhs, "nclust",
      "otherwise uniform on 1..ntime-1 occasions")
  for (k in seq_along(matched)) {
    s   <- st$draws[[matched[k]]]
    tgt <- s$target
    can <- sub("_raw$", "", sub("^d\\$", "", tgt))
    nm  <- if (s$fn == "sample") "size" else "n"
    nval <- as.numeric(dyn[[k]]$args[[nm]])
    lvl <- if (nval == N && n_obs != N) "subject" else if (nval == n_obs) "observation" else "?"
    if (tgt == "nclust") lvl <- "subject"
    std_rhs <- vapply(Filter(function(z) z$target %in% c(tgt, can, paste0("d$", can)),
                             st$assigns), `[[`, "", "rhs")
    std <- any(grepl("scale(", std_rhs, fixed = TRUE))
    add("design", sprintf("draw_%02d", k), s$call, can,
        sprintf("level=%s; standardised=%s; assigned_as=%s", lvl, std, tgt))
  }
  used <- unique(unlist(lapply(c(list(mf), if (dgp == "GenBinaryY") list(a$lv.formula, a$t.formula)
                                 else list(a$re.formula)), all.vars)))
  cols <- setdiff(names(dat), c("id", "time", "binY", "sim_id"))
  add("design", "data_columns", paste(names(dat), collapse = ", "))
  unused <- setdiff(cols, used)
  add("design", "unused_covariates", if (length(unused)) paste(unused, collapse = ", ") else "none",
      "", "drawn (consume random numbers) but enter no part of the DGM")
  add("rng", "seed_rule", seed_txt, "", "set.seed() per replicate in run_full.R")
  add("rng", "RNGkind", WORKER_RNGKIND, "",
      sprintf("set by crew %s / mirai %s in every worker; NOT set by run_full.R (plain R session: %s)",
              packageVersion("crew"), packageVersion("mirai"), SESSION_RNGKIND))
  do.call(rbind, rows)
}

all_rows <- list(); invariance <- list()
for (dgp in DGPS) for (model in 1:5) for (arm in names(ARMS)) {
  caps <- lapply(seq_len(nrow(CELLS)), function(r)
    capture_one(dgp, model, arm, CELLS$N[r], CELLS$ntime[r], CELLS$i[r]))
  key  <- function(x) x[x$block != "derived", c("parameter", "value", "terms", "note")]
  same <- all(vapply(caps[-1], function(x) identical(key(x), key(caps[[1]])), logical(1)))
  invariance[[length(invariance) + 1L]] <- data.frame(dgp = dgp, model = model, arm = arm,
                                                      invariant_over_cells = same)
  all_rows[[length(all_rows) + 1L]] <- caps[[1]]
}
tab <- do.call(rbind, all_rows)
inv <- do.call(rbind, invariance)
write.csv(tab, OUT, row.names = FALSE)
write.csv(inv, sub("[.]csv$", "_invariance.csv", OUT), row.names = FALSE)
cat(sprintf("rundir: %s\nrows: %d | generator x arm combinations: %d | invariant over 3 cells: %d/%d\n",
            RUNDIR, nrow(tab), nrow(unique(tab[, c("dgp", "model", "arm")])),
            sum(inv$invariant_over_cells), nrow(inv)))
cat("written:", OUT, "\n")
