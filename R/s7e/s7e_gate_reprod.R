# s7e_gate_reprod.R -- GF4 (a): three chunks of the S7e run regenerated in a new process must be identical to the run's files
# (records, study rows, convergence-check fits, GEE rows; timings excluded).
source(path.expand("~/MMM/R/s7e/s7e_common.R"))
TK <- s7e_tasks()
pick <- TK[(scen == "M5" & N == 145L & ntime == 10L & label == "manual" & chunk == 3L) | (scen == "M1" & N == 10L & ntime == 2L & label == "mmlb" & chunk == 1L) |
           (scen == "M3" & N == 55L & ntime == 4L & label == "manual" & chunk == 10L)]
stopifnot(nrow(pick) == 3L)
res <- rbindlist(lapply(seq_len(nrow(pick)), function(k) {
  t <- pick[k]; old <- readRDS(t$file); new <- run_chunk_s7e(t, out_file = NULL)
  ok <- vapply(c("recs", "study", "cc", "gee"), function(x) isTRUE(all.equal(old[[x]], new[[x]], tolerance = 0, check.attributes = FALSE)), TRUE)
  data.table(chunk = basename(t$file), recs = ok[["recs"]], study = ok[["study"]], cc = ok[["cc"]], gee = ok[["gee"]])
}))
res[, pass := recs & study & cc & gee]; print(res); fwrite(res, file.path(TABDIR, "V4_S7e_GF4a.csv"))
