# s7e_run.R [max_chunks] -- S7e production run, resumable (existing chunk files are skipped). Worker count read at every
# scheduling step from results_v4/S7e/WORKERS (capped at 28, the thermal ceiling of the host; guard R/s7e/s7e_thermal_guard.sh).
# Largest chunks first (M5, then by N * T). Optional argument: run at most this many chunks (dry run), in the order given by DRY_KEYS if set.
source(path.expand("~/MMM/R/s7e/s7e_common.R"))
a <- commandArgs(TRUE); MAXC <- if (length(a)) as.integer(a[1]) else NA_integer_
CTL <- file.path(OUT_S7E, "WORKERS"); LOG <- file.path(OUT_S7E, "run_progress.log")
TK <- s7e_tasks(); TK <- TK[!file.exists(file)][order(-(model == 5L), -N * ntime)]
if (!is.na(MAXC)) { dk <- Sys.getenv("DRY_KEYS"); if (nzchar(dk)) TK <- TK[basename(file) %in% strsplit(dk, ",")[[1]]]; TK <- head(TK, MAXC) }
cat(format(Sys.time()), "chunks to run:", nrow(TK), "\n", file = LOG, append = TRUE)
nw <- function() { w <- suppressWarnings(as.integer(readLines(CTL, n = 1L))); if (length(w) != 1L || is.na(w)) 1L else max(0L, min(w, 28L)) }
active <- list(); q <- seq_len(nrow(TK)); ok <- 0L; bad <- 0L; t0 <- Sys.time(); last <- t0
while (length(q) || length(active)) {
  W <- nw()
  while (length(active) < W && length(q)) {
    r <- q[1]; q <- q[-1]
    p <- parallel::mcparallel({ data.table::setDTthreads(1L)
      tryCatch({ run_chunk_s7e(TK[r]); 0L }, error = function(e) { message("TASK FAILED ", basename(TK$file[r]), ": ", conditionMessage(e)); 1L }) })
    active[[as.character(p$pid)]] <- p
  }
  if (length(active)) {
    res <- parallel::mccollect(active, wait = FALSE, timeout = 5)
    if (length(res)) for (k in names(res)) { if (identical(res[[k]], 0L)) ok <- ok + 1L else bad <- bad + 1L; active[[k]] <- NULL }
  } else Sys.sleep(5)
  if (difftime(Sys.time(), last, units = "secs") > 300) { last <- Sys.time()
    cat(sprintf("%s workers %d | done %d ok %d failed | left %d | %.2f h elapsed\n", format(Sys.time()), length(active), ok, bad, length(q),
                as.numeric(difftime(Sys.time(), t0, units = "hours"))), file = LOG, append = TRUE) }
}
cat(sprintf("%s %s: %d ok, %d failed; %.2f h\n", format(Sys.time()), if (is.na(MAXC)) "DONE" else "DRY FINISHED", ok, bad,
            as.numeric(difftime(Sys.time(), t0, units = "hours"))), file = LOG, append = TRUE)
