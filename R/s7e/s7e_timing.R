# s7e_timing.R -- per-replicate cost of the full S7e pipeline on copula data (no performance measure is computed or printed).
source(path.expand("~/MMM/R/s7e/s7e_common.R"))
cells <- data.table(scen = c("M1", "M1", "M3", "M5", "M5", "M2", "M4"), N = c(10L, 145L, 85L, 145L, 25L, 55L, 145L), ntime = c(2L, 10L, 6L, 10L, 3L, 4L, 8L))
cells[, `:=`(label = "manual", chunk = 1L, model = as.integer(sub("M", "", scen)), file = NA_character_)]
res <- rbindlist(lapply(seq_len(nrow(cells)), function(k) {
  t <- cells[k]; t0 <- Sys.time()
  r <- run_chunk_s7e(t, ids = 1:8, out_file = NULL)
  data.table(t[, .(scen, N, ntime)], secs_per_rep = as.numeric(difftime(Sys.time(), t0, units = "secs")) / 8,
             t(colMeans(r$secs)), n_study_rows = nrow(r$study), n_gee_rows = nrow(r$gee), n_cc = if (is.null(r$cc)) 0L else nrow(r$cc))
}))
print(res, digits = 3)
fwrite(res, file.path(OUT_S7E, "timing_probe.csv"))
