# =============================================================================
# s10_degen_run.R [max_chunks] -- S10 "Degeneri": resumable production run of s10_degen_common.R::dg_chunk over the
# 47,520 chunks (existing output files skipped). Worker count read at every scheduling step from
# results_v4/S10/degen/WORKERS (capped at 14). Thermal guard inside the scheduler: every 240 s Tctl is read from
# `sensors`; Tctl > 90 C -> WORKERS - 2 (min 4), Tctl < 82 C -> WORKERS + 1 (max 14); running chunks are never killed.
# A file results_v4/S10/degen/STOP stops scheduling (running chunks finish). Logs in results_v4/S10/degen/logs/.
# =============================================================================
source(path.expand("~/MMM/R/s10/degen/s10_degen_common.R"))
a <- commandArgs(TRUE); MAXC <- if (length(a)) as.integer(a[1]) else NA_integer_
CTL <- file.path(OUT_DG, "WORKERS"); STOPF <- file.path(OUT_DG, "STOP")
LOG <- file.path(OUT_DG, "logs/run_progress.log"); TLOG <- file.path(OUT_DG, "logs/thermal.log")
WMAX <- 14L; WMIN <- 4L
if (!file.exists(CTL)) writeLines("10", CTL)
TK <- dg_tasks(); TK <- TK[!file.exists(out)][order(-N * ntime, src, dgp, scen, N, ntime, label, chunk)]
if (!is.na(MAXC)) { sel <- Sys.getenv("DG_SEL"); if (nzchar(sel)) TK <- TK[eval(parse(text = sel))]; TK <- head(TK, MAXC) }
for (g in unique(TK$dgp)) invisible(dg_env(g))   # environments loaded once in the parent, inherited by the forks
tctl <- function() { x <- tryCatch(system("sensors 2>/dev/null", intern = TRUE), error = function(e) character(0)); x <- grep("Tctl", x, value = TRUE)
  if (!length(x)) NA_real_ else as.numeric(sub("^.*\\+([0-9.]+).*$", "\\1", x[1])) }
nw <- function() { w <- suppressWarnings(as.integer(readLines(CTL, n = 1L))); if (length(w) != 1L || is.na(w)) WMIN else max(0L, min(w, WMAX)) }
lavg <- function() tryCatch(strsplit(readLines("/proc/loadavg"), " ")[[1]][1], error = function(e) NA)
cat(format(Sys.time()), "chunks to run:", nrow(TK), "| Tctl", tctl(), "| load", lavg(), "\n", file = LOG, append = TRUE)
active <- list(); q <- seq_len(nrow(TK)); ok <- 0L; bad <- 0L; t0 <- Sys.time(); last <- t0 - 1000; lastT <- t0 - 1000
while (length(q) || length(active)) {
  if (file.exists(STOPF) && length(q)) { cat(format(Sys.time()), "STOP file found: no new chunks\n", file = LOG, append = TRUE); q <- integer(0) }
  if (difftime(Sys.time(), lastT, units = "secs") > 240) { lastT <- Sys.time(); T <- tctl(); W <- nw(); Wn <- W
    if (is.finite(T) && T > 90 && W > WMIN) Wn <- max(WMIN, W - 2L) else if (is.finite(T) && T < 82 && W < WMAX) Wn <- W + 1L
    if (Wn != W) writeLines(as.character(Wn), CTL)
    cat(sprintf("%s Tctl %.1f | load %s | WORKERS %d -> %d | active %d\n", format(Sys.time()), T, lavg(), W, Wn, length(active)), file = TLOG, append = TRUE) }
  W <- nw()
  while (length(active) < W && length(q)) {
    r <- q[1]; q <- q[-1]
    p <- parallel::mcparallel({ data.table::setDTthreads(1L)
      tryCatch({ dg_chunk(TK[r]); 0L }, error = function(e) { message("TASK FAILED ", basename(TK$out[r]), ": ", conditionMessage(e)); 1L }) })
    active[[as.character(p$pid)]] <- p
  }
  if (length(active)) {
    res <- parallel::mccollect(active, wait = FALSE, timeout = 5)
    if (length(res)) for (k in names(res)) { if (identical(res[[k]], 0L)) ok <- ok + 1L else bad <- bad + 1L; active[[k]] <- NULL }
  } else Sys.sleep(5)
  if (difftime(Sys.time(), last, units = "secs") > 300) { last <- Sys.time()
    cat(sprintf("%s workers %d | done %d ok %d failed | left %d | %.2f h elapsed | Tctl %.1f | load %s\n", format(Sys.time()), length(active), ok, bad, length(q),
                as.numeric(difftime(Sys.time(), t0, units = "hours")), tctl(), lavg()), file = LOG, append = TRUE) }
}
cat(sprintf("%s %s: %d ok, %d failed; %.2f h\n", format(Sys.time()), if (is.na(MAXC)) "DONE" else "DRY FINISHED", ok, bad,
            as.numeric(difftime(Sys.time(), t0, units = "hours"))), file = LOG, append = TRUE)
