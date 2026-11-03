# s7g_smoke.R -- S7g tool check before pre-registration: 2 replicates on 5 design points, nothing saved; timings per fit.
source(path.expand("~/MMM/R/s7g/s7g_common.R"))
TK <- s7g_tasks()
pick <- TK[chunk == 1L & ((scen == "M2" & N == 145L & ntime == 10L) | (scen == "M1c" & N == 45L & ntime == 6L & dgp == "GenBinaryY"))]
res <- rbindlist(lapply(seq_len(nrow(pick)), function(k) { t <- pick[k]; t0 <- Sys.time()
  r <- run_chunk_s7g(t, ids = c(1L, 101L), out_file = NULL)$rows; r[, wall_task := as.numeric(difftime(Sys.time(), t0, units = "secs"))]; r }))
print(res[, .(dgp, scen, N, ntime, sim_id, nAGQ, estimate = signif(estimate, 6), se = signif(se, 6), D = signif(D, 4), status = substr(status, 1, 40), secs)])
print(res[, .(mean_secs = mean(secs)), by = .(nAGQ, big = N * ntime >= 1000)])
