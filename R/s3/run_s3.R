# =============================================================================
# run_s3.R -- S3 production run (18,720 chunks), resumable: chunks whose RDS
# exists are skipped. Same orchestration as R/run_full.R (crew, 50 workers,
# R_LIBS and one BLAS thread exported before the controller starts, bounded
# backlog). Log: results_v4/S3/run.log (stable path).
# Env: MMM_S3_BLOCKS (e.g. "F1,F2"), MMM_S3_KEYS (file with keys, pilot).
# =============================================================================
source(path.expand("~/MMM/R/s3/s3_common.R"))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep),
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
suppressMessages(library(crew))
dir.create(OUT_S3, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT_S3, "run.log")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(format(Sys.time(), "[%F %T] "), m, "\n", file = LOG, append = TRUE, sep = "") }

tasks <- s3_tasks()
if (nzchar(Sys.getenv("MMM_S3_BLOCKS"))) tasks <- tasks[tasks$block %in% strsplit(Sys.getenv("MMM_S3_BLOCKS"), ",")[[1]], ]
if (nzchar(Sys.getenv("MMM_S3_KEYS"))) tasks <- tasks[tasks$key %in% readLines(Sys.getenv("MMM_S3_KEYS")), ]
tasks <- tasks[order(-tasks$N * tasks$ntime, tasks$block, tasks$key), ]   # longest first
done <- file.exists(tasks$file); todo <- tasks[!done, ]
say("S3 run: total %d | done %d | to run %d | commit %s", nrow(tasks), sum(done), nrow(todo),
    system("git -C ~/MMM log --oneline -n 1", intern = TRUE))
if (!nrow(todo)) quit(status = 0)

worker <- function(t, repo) { source(file.path(repo, "R/s3/s3_common.R")); run_chunk_s3(t) }
WORKERS <- as.integer(Sys.getenv("MMM_S3_WORKERS", "50"))
ctrl <- crew_controller_local(workers = WORKERS, seconds_idle = 60, tasks_max = 50, garbage_collection = TRUE)
ctrl$start()
push_one <- function(k) {
  t <- as.list(todo[k, ])
  ctrl$push(command = FUN(t, repo), data = list(FUN = worker, t = t, repo = REPO), name = t$key)
}
BACKLOG <- 400L; pushed <- 0L
while (pushed < min(BACKLOG, nrow(todo))) { pushed <- pushed + 1L; push_one(pushed) }
t0 <- Sys.time(); ndone <- 0L; nerr <- 0L; last <- Sys.time()
while (ndone + nerr < nrow(todo)) {
  ctrl$wait(mode = "one", seconds_timeout = 600)
  repeat {
    p <- ctrl$pop(); if (is.null(p)) break
    e <- p$error[[1]]
    if (!is.null(e) && !is.na(e)) { nerr <- nerr + 1L; say("ERROR %s :: %s", p$name, substr(as.character(e), 1, 160)) }
    else { ndone <- ndone + 1L; cat(sprintf("%s\t%.1f\n", p$name, p$seconds), file = file.path(OUT_S3, "chunk_seconds.tsv"), append = TRUE) }
    if (pushed < nrow(todo)) { pushed <- pushed + 1L; push_one(pushed) }
  }
  if (as.numeric(difftime(Sys.time(), last, units = "secs")) > 300) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); tot <- ndone + nerr
    say("%d/%d chunks (%.1f%%) | %d errors | %.1f min | ETA %.0f min", tot, nrow(todo), 100 * tot / nrow(todo),
        nerr, el, if (tot > 0) el * (nrow(todo) - tot) / tot else NA)
    last <- Sys.time()
  }
}
ctrl$terminate()
say("DONE: %d ok, %d errors, %.1f min", ndone, nerr, as.numeric(difftime(Sys.time(), t0, units = "mins")))
