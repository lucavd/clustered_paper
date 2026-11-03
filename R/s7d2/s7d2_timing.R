# s7d2_timing.R -- time pilot for organising the S7d2 run (no measure of this script enters the paper or the predictions):
# wall-clock of one_rep_s7d2() on replicates 1-10 (H1) at 36 points (GenBinaryY M1, GLMM-AR1 M5, copula M3; N 10, 25, 65,
# 145; T 2, 6, 10), one point per worker; prints only times and the number of test evaluations.
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))
P <- CJ(k = 1:3, N = c(10L, 25L, 65L, 145L), ntime = c(2L, 6L, 10L))
P[, `:=`(dgp = c("GenBinaryY", "GLMM-AR1", DGP_COP)[k], scen = c("M1", "M5", "M3")[k], model = c(1L, 5L, 3L)[k], label = "mmlb")]
one <- function(r) { t <- P[r]; gen <- s7d2_gen(t); th <- truth_of(t$dgp, t$scen, t$ntime, t$label)
  t0 <- Sys.time(); x <- rbindlist(lapply(1:10, function(i) one_rep_s7d2(gen, t, i, th)), fill = TRUE)
  data.table(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, wall_per_rep = as.numeric(difftime(Sys.time(), t0, units = "secs")) / 10,
             ci_secs = mean(x$secs, na.rm = TRUE), n_eval = mean(x$n_eval, na.rm = TRUE), n_eval_max = max(x$n_eval, na.rm = TRUE), secs_max = max(x$secs, na.rm = TRUE), n_err = sum(!is.na(x$ci_error))) }
R <- rbindlist(parallel::mclapply(seq_len(nrow(P)), one, mc.cores = 18, mc.preschedule = FALSE))
print(R, digits = 3)
fwrite(R, file.path(OUT_S7D2, "timing_pilot.csv"))
