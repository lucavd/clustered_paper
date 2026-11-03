# s6b_ns_run.R -- S6b part B production run (D-34): 4 scenarios x 16 N x 6 T x 2 arms x 1,000 replicates, in chunks
# of 100 (7,680 files under results_v4/S6b/ns/rows/). Resumable: existing chunk files are skipped.
# The number of concurrent workers is read at every scheduling step from results_v4/S6b/ns/WORKERS (one integer),
# so the run can share the host with other jobs under the 50-core cap and be widened without a restart.
# Usage: Rscript --vanilla R/s6b/s6b_ns_run.R
source(path.expand("~/MMM/R/s6b/s6b_ns_common.R"))
CTL <- file.path(OUT_NS, "WORKERS"); LOG <- file.path(OUT_NS, "run_progress.log")
TK <- ns_tasks(); TK <- TK[!file.exists(file)][order(-N * ntime)]
cat(format(Sys.time()), "chunks to run:", nrow(TK), "\n", file = LOG, append = TRUE)
nw <- function() { w <- suppressWarnings(as.integer(readLines(CTL, n = 1L))); if (length(w) != 1L || is.na(w)) 1L else max(0L, min(w, 50L)) }
active <- list(); q <- seq_len(nrow(TK)); ok <- 0L; bad <- 0L; t0 <- Sys.time(); last <- t0
while (length(q) || length(active)) {
  W <- nw()
  while (length(active) < W && length(q)) {
    r <- q[1]; q <- q[-1]
    p <- parallel::mcparallel({ data.table::setDTthreads(1L)
      tryCatch({ run_chunk_ns(TK[r]); 0L }, error = function(e) { message("TASK FAILED ", basename(TK$file[r]), ": ", conditionMessage(e)); 1L }) })
    active[[as.character(p$pid)]] <- p
  }
  if (length(active)) {
    res <- parallel::mccollect(active, wait = FALSE, timeout = 5)
    if (length(res)) { for (k in names(res)) { if (identical(res[[k]], 0L)) ok <- ok + 1L else bad <- bad + 1L; active[[k]] <- NULL } }
  } else Sys.sleep(5)
  if (difftime(Sys.time(), last, units = "secs") > 300) { last <- Sys.time()
    el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    cat(sprintf("%s workers %d | done %d ok %d failed | left %d | %.2f h elapsed\n", format(Sys.time()), length(active), ok, bad, length(q), el), file = LOG, append = TRUE) }
}
cat(sprintf("%s DONE: %d ok, %d failed; %.2f h\n", format(Sys.time()), ok, bad, as.numeric(difftime(Sys.time(), t0, units = "hours"))), file = LOG, append = TRUE)
