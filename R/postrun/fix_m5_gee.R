
# =============================================================================
# fix_m5_gee.R -- recompute ONLY the uncorrected GEE under scenario 5.
#
# Finding A1: do_analysis_gee() routed model 5 to analyze_model4_gee(), which
# fits binY ~ time*trt + age -- omitting sex and steroids. The MMM and the GLMM
# were both repaired on 2026-08-05; the GEE arm was not, and analyze_model5_gee()
# has sat there unused. Note gee_ss was NOT affected: it refits from `mean_rhs`,
# which is correct for M5, which is exactly why `gee` and `gee_md` disagree at
# model 5 and agree everywhere else (verified: 0/1000 vs 1000/1000 identical).
#
# Scope: technique "gee", model 5 only. 3 DGP x 16 N x 6 T x 2 arms x 10 chunks.
# The design is a deterministic function of the seed, and the seed is intact, so
# the SAME datasets are regenerated -- the other techniques' stored rows remain
# valid and are untouched.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
suppressMessages({library(crew); library(dplyr)})

SRC    <- path.expand("~/MMM/R")
OUTDIR <- path.expand("~/MMM/results_m5gee")
dir.create(OUTDIR, showWarnings = FALSE)

N_LEVELS <- c(10, 15, 20, seq(25, 145, 10)); T_LEVELS <- c(2, 3, 4, 6, 8, 10)
LABELS <- c("mmlb", "manual"); DGPS <- c("GenBinaryY", "GLMM-native", "GLMM-AR1")
N_SIM <- 1000; CHUNK <- 100; MODEL <- 5L

cells <- expand.grid(dgp = DGPS, N = N_LEVELS, ntime = T_LEVELS, label = LABELS,
                     chunk = seq_len(N_SIM / CHUNK), stringsAsFactors = FALSE)
cells$key  <- sprintf("%s_m5_N%d_T%d_%s_c%02d", gsub("[^A-Za-z0-9]", "", cells$dgp),
                      cells$N, cells$ntime, cells$label, cells$chunk)
cells$file <- file.path(OUTDIR, paste0(cells$key, ".rds"))
cells <- cells[order(cells$dgp, cells$N, cells$ntime, cells$label, cells$chunk), ]
todo <- cells[!file.exists(cells$file), ]
cat(sprintf("chunks: %d | done: %d | to run: %d\n", nrow(cells),
            sum(file.exists(cells$file)), nrow(todo)))
if (!nrow(todo)) { cat("nothing to do\n"); quit(status = 0) }

run_chunk <- function(dgp, N, ntime, label, chunk, src, lib, chunk_size, out_file) {
  .libPaths(lib)
  suppressMessages({library(geepack); library(dplyr); library(tibble); library(MASS)
                    library(MMLB); library(fastGHQuad)})
  sp <- function(f) sys.source(file.path(src, f), envir = globalenv())
  for (f in c("analysis_model5.R")) sp(f)
  msm <- file.path(dirname(src), "MMM_sim_manual")
  if (dgp == "GenBinaryY") {
    for (f in c("generation_model1.R","generation_model2.R","generation_model3.R",
                "generation_model4.R","generation_model5.R","generate_data.R",
                "make_all_data.R")) sp(f)
  } else if (dgp == "GLMM-native") {
    sys.source(file.path(msm, "gen_glmm_data.R"), envir = globalenv())
    for (f in c("generation_model1.R","generation_model2.R","generation_model3.R",
                "generation_model4.R","generation_model5.R","generate_data.R",
                "make_all_data.R"))
      sys.source(file.path(msm, "R_manual", f), envir = globalenv())
  } else {
    sys.source(file.path(msm, "gen_ar1_data.R"), envir = globalenv())
    sys.source(file.path(msm, "R_ar1", "generation_models_ar1.R"), envir = globalenv())
    sys.source(file.path(msm, "R_manual", "generate_data.R"), envir = globalenv())
    sys.source(file.path(msm, "R_ar1", "make_all_data.R"), envir = globalenv())
  }
  base <- (chunk - 1L) * chunk_size
  acc <- vector("list", chunk_size)
  for (j in seq_len(chunk_size)) {
    i <- base + j
    set.seed(20260805L + 1000L * 5L + i)          # identical to the production seed
    d <- tryCatch(make_all_data(N, ntime, i, 5L, label), error = function(e) NULL)
    if (is.null(d)) next
    out <- tryCatch(analyze_model5_gee(d, N, ntime, label), error = function(e) NULL)
    if (is.null(out)) {
      out <- data.frame(term = NA_character_, estimate = NA_real_, se = NA_real_,
                        stat = NA_real_, p = NA_real_, method = "gee",
                        N = N, ntime = ntime, label = label,
                        converged = FALSE, status = "fit error")
    } else { out$converged <- TRUE; out$status <- NA_character_ }
    out$model <- 5L; out$sim_id <- i; out$dgp <- dgp
    acc[[j]] <- out
  }
  res <- dplyr::bind_rows(acc)
  tmp <- paste0(out_file, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out_file)
  nrow(res)
}

ctrl <- crew_controller_local(workers = 50, seconds_idle = 30)
ctrl$start()
for (k in seq_len(nrow(todo))) {
  ctrl$push(command = run_chunk(dgp, N, ntime, label, chunk, src, lib, chunk_size, out_file),
            data = list(dgp = todo$dgp[k], N = todo$N[k], ntime = todo$ntime[k],
                        label = todo$label[k], chunk = todo$chunk[k], src = SRC,
                        lib = LIB, chunk_size = CHUNK, out_file = todo$file[k],
                        run_chunk = run_chunk),
            name = todo$key[k])
}
cat("pushed", nrow(todo), "chunks\n"); flush.console()
t0 <- Sys.time(); nerr <- 0L
while (!ctrl$empty()) {
  r <- ctrl$pop(scale = TRUE)
  if (is.null(r)) { Sys.sleep(0.4); next }
  if (!is.null(r$error) && !is.na(r$error)) { nerr <- nerr + 1L
    if (nerr <= 5) cat("ERROR:", r$name, "|", substr(r$error, 1, 160), "\n") }
  n <- length(list.files(OUTDIR, pattern = "[.]rds$"))
  if (n %% 100 == 0) { el <- as.numeric(difftime(Sys.time(), t0, units = "min"))
    cat(sprintf("[%s] %d/%d  %.1f min\n", format(Sys.time(), "%H:%M"), n,
                nrow(cells), el)); flush.console() }
}
ctrl$terminate()
cat(sprintf("DONE files=%d errors=%d wall=%.1f min\n",
            length(list.files(OUTDIR, pattern = "[.]rds$")), nerr,
            as.numeric(difftime(Sys.time(), t0, units = "min"))))
