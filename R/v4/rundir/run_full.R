# =============================================================================
# run_full.R -- full simulation rerun, three DGPs, extended grid.
# Launched 2026-08-05 after the pre-launch validation suite passed.
#
# Design notes:
#  * crew with 100 workers (user instruction; the box is shared, 128 cores).
#  * R_LIBS is exported before the controller starts: crew workers are fresh R
#    processes that do NOT inherit .libPaths(), and without this they cannot
#    load MMLB/lme4 and the controller hangs with no error (diagnosed today).
#  * Work is chunked into (cell x 100 replicates) tasks and each completed
#    chunk is written to its own RDS immediately. The run is therefore
#    RESUMABLE: relaunching skips chunks whose RDS already exists. The host
#    dropped off the network mid-job earlier today, so this is not theoretical.
#  * Each worker sources its own copy of the three DGP families into separate
#    environments, so the three generators never collide by name.
# =============================================================================

LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

# One BLAS thread per worker. R links OpenBLAS-pthread here, which by default
# spawns a thread per core INSIDE each worker: with 50 workers that measured 95
# threads each and a load average of 287 on 128 cores -- 5.7x the agreed cap on
# a shared machine. The work is embarrassingly parallel ACROSS replicates, so
# intra-fit threading buys nothing and only oversubscribes. Workers are fresh R
# processes and do not inherit these unless they are set in the environment
# before the controller starts, exactly as with R_LIBS above.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           RCPP_PARALLEL_NUM_THREADS = "1")

suppressMessages({library(crew); library(dplyr)})

SRC     <- normalizePath(".")
# Results live OUTSIDE the per-job scratch directory: each dispatch gets a new
# workdir, so keeping them in "./results" would orphan completed chunks on every
# relaunch. This path persists across jobs and across reboots.
OUTDIR  <- path.expand("~/MMM/results_full")
dir.create(OUTDIR, showWarnings = FALSE)

N_LEVELS <- c(10, 15, 20, seq(25, 145, 10))
T_LEVELS <- c(2, 3, 4, 6, 8, 10)
# MMM_MODELS restricts the run to a subset of scenarios, so a redesign of
# some scenarios can be regenerated without touching the others.
MODELS   <- if (nzchar(Sys.getenv("MMM_MODELS")))
  as.integer(strsplit(Sys.getenv("MMM_MODELS"), ",")[[1]]) else 1:5
LABELS   <- c("mmlb", "manual")
DGPS     <- c("GenBinaryY", "GLMM-native", "GLMM-AR1")
N_SIM    <- 1000
CHUNK    <- 100

# Smoke mode: exercise the full orchestration path (crew, all three DGPs,
# incremental writes, resumability) on a few cells before committing hours.
if (nzchar(Sys.getenv("MMM_SMOKE"))) {
  N_LEVELS <- c(10, 45); T_LEVELS <- c(2, 6); MODELS <- c(1, 5)
  N_SIM <- 200; CHUNK <- 100
  OUTDIR <- file.path(SRC, "results_smoke"); dir.create(OUTDIR, showWarnings = FALSE)
  cat("*** SMOKE MODE ***\n")
}
TECHS    <- c("mmm", "gee", "glmer", "gee_ss", "glmer_lrt")

cells <- expand.grid(dgp = DGPS, model = MODELS, N = N_LEVELS, ntime = T_LEVELS,
                     label = LABELS, chunk = seq_len(N_SIM / CHUNK),
                     stringsAsFactors = FALSE)
cells$key <- sprintf("%s_m%d_N%d_T%d_%s_c%02d",
                     gsub("[^A-Za-z0-9]", "", cells$dgp), cells$model,
                     cells$N, cells$ntime, cells$label, cells$chunk)
cells$file <- file.path(OUTDIR, paste0(cells$key, ".rds"))

# Order tasks so a cell's 10 chunks run CONSECUTIVELY. expand.grid varies its
# first column fastest and `chunk` slowest, which meant the queue swept all 2880
# cells at chunk 1 before returning for chunk 2 -- so no cell reached its full
# 1000 replicates until ~90% of the run was done, and partial tables showed
# 100-replicate cells throughout. Sorting by cell makes completed cells
# accumulate steadily, so the tables are readable while the run proceeds.
cells <- cells[order(cells$dgp, cells$model, cells$N, cells$ntime,
                     cells$label, cells$chunk), ]

done <- file.exists(cells$file)
cat(sprintf("total chunks: %d | already done: %d | to run: %d\n",
            nrow(cells), sum(done), sum(!done)))
todo <- cells[!done, ]
if (!nrow(todo)) { cat("nothing to do\n"); quit(status = 0) }

# ---- the unit of work, evaluated in a worker -------------------------------
run_chunk <- function(dgp, model, N, ntime, label, chunk, src, lib, techs,
                      chunk_size, out_file) {
  .libPaths(lib)
  suppressMessages({
    library(MMLB); library(lme4); library(geepack); library(geesmv)
    library(dplyr); library(tibble); library(MASS); library(fastGHQuad)
  })
  sp <- function(f) sys.source(file.path(src, f), envir = globalenv())
  for (f in c("analysis_model1.R", "analysis_model2.R", "analysis_model3.R",
              "analysis_model4.R", "analysis_model5.R", "analysis_smallsample.R",
              "do_analysis.R")) sp(f)
  if (dgp == "GenBinaryY") {
    for (f in c("generation_model1.R", "generation_model2.R", "generation_model3.R",
                "generation_model4.R", "generation_model5.R",
                "generate_data.R", "make_all_data.R")) sp(f)
  } else if (dgp == "GLMM-native") {
    sp("gen_glmm_data.R")
    for (f in c("manual/generation_model1.R", "manual/generation_model2.R",
                "manual/generation_model3.R", "manual/generation_model4.R",
                "manual/generation_model5.R", "manual/generate_data.R",
                "manual/make_all_data.R")) sp(f)
  } else {
    sp("gen_ar1_data.R")
    for (f in c("ar1/generation_models_ar1.R", "manual/generate_data.R",
                "ar1/make_all_data.R")) sp(f)
  }

  base <- (chunk - 1L) * chunk_size
  acc <- vector("list", chunk_size)
  for (j in seq_len(chunk_size)) {
    i <- base + j
    set.seed(20260805L + 1000L * model + i)
    d <- tryCatch(make_all_data(N, ntime, i, model, label), error = function(e) NULL)
    if (is.null(d)) next
    rows <- lapply(techs, function(tc)
      tryCatch(do_analysis(d, N, ntime, label, model, tc, i), error = function(e) NULL))
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) acc[[j]] <- dplyr::bind_rows(rows)
  }
  res <- dplyr::bind_rows(Filter(Negate(is.null), acc))
  if (nrow(res)) {
    res$dgp <- dgp
    tmp <- paste0(out_file, ".tmp")
    saveRDS(res, tmp, compress = "xz")
    file.rename(tmp, out_file)   # atomic: a reader never sees a partial file
  }
  c(rows = nrow(res))
}

# 50 workers (user instruction 2026-08-05, after the box became unresponsive at
# 100). The real concurrency is HIGHER than the worker count: analyze_gee_smallsample
# forks a child per correction for its hard timeout, so 100 workers drove load
# average to ~161 on 128 cores. At 50 the peak stays under the core count even
# with forks outstanding. This machine also runs an interactive desktop session,
# so leaving headroom is what keeps it responsive.
ctrl <- crew_controller_local(workers = 50, seconds_idle = 60,
                              tasks_max = 50, garbage_collection = TRUE)
ctrl$start()
cat("controller started, pushing", nrow(todo), "tasks\n")

# Push in bounded waves rather than queueing all 28k tasks at once: the launcher
# held ~660 MB of queued task payloads in the previous attempt, and a huge queue
# makes an interrupted run harder to reason about.
push_one <- function(k) {
  r <- todo[k, ]
  # NB: crew evaluates `command` in a FRESH worker session that does not inherit
  # this script's globals -- pushing `run_chunk(...)` alone fails with
  # "could not find function run_chunk" (caught by the smoke test). The function
  # itself must travel in `data`.
  ctrl$push(
    command = FUN(dgp, model, N, ntime, label, chunk, src, lib, techs,
                  chunk_size, out_file),
    data = list(FUN = run_chunk,
                dgp = r$dgp, model = r$model, N = r$N, ntime = r$ntime,
                label = r$label, chunk = r$chunk, src = SRC, lib = LIB,
                techs = TECHS, chunk_size = CHUNK, out_file = r$file),
    name = r$key
  )
}
BACKLOG <- 400L
next_k  <- 0L
pushed  <- 0L
while (pushed < min(BACKLOG, nrow(todo))) { pushed <- pushed + 1L; push_one(pushed) }
next_k <- pushed
cat(sprintf("primed %d tasks (backlog cap %d) at %s\n", pushed, BACKLOG, format(Sys.time())))

t0 <- Sys.time(); ndone <- 0L; nerr <- 0L; last <- Sys.time()
while (ndone + nerr < nrow(todo)) {
  ctrl$wait(mode = "one", seconds_timeout = 600)
  repeat {
    p <- ctrl$pop()
    if (is.null(p)) break
    e <- p$error[[1]]
    if (!is.null(e) && !is.na(e)) {
      nerr <- nerr + 1L
      cat(sprintf("ERROR %s :: %s\n", p$name, substr(as.character(e), 1, 160)))
    } else ndone <- ndone + 1L
    if (next_k < nrow(todo)) { next_k <- next_k + 1L; push_one(next_k) }
  }
  if (as.numeric(difftime(Sys.time(), last, units = "secs")) > 120) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    tot <- ndone + nerr
    cat(sprintf("[%s] %d/%d chunks (%.1f%%) | %d errors | %.1f min elapsed | ETA %.0f min\n",
                format(Sys.time(), "%H:%M:%S"), tot, nrow(todo),
                100 * tot / nrow(todo), nerr, el,
                if (tot > 0) el * (nrow(todo) - tot) / tot else NA))
    flush.console(); last <- Sys.time()
  }
}
ctrl$terminate()
cat(sprintf("\nDONE: %d chunks ok, %d errors, %.1f min\n", ndone, nerr,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("result files:", length(list.files(OUTDIR, pattern = "[.]rds$")), "\n")
