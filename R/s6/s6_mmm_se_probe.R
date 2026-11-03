# =============================================================================
# s6_mmm_se_probe.R -- S6 follow-up of the PI's question on the marginalized model (MMM).
# No new simulated data: the first 200 replicates of three existing cells are regenerated with
# their own seeds (set.seed(20260805 + 1000 * model + i), L'Ecuyer-CMRG, as R/v4/rundir/run_full.R)
# and refitted with mm() exactly as in R/v4/rundir/analysis_model*.R. For time:trt the probe keeps
# the model-based SE (mod.cov, the one the study uses, via summary()) and the robust SE that mm()
# also returns (rob.cov), and checks that the estimate equals the stored one (results_v4/S4/rows_v4.rds).
# Cells: GenBinaryY, H0, N = 145, T = 10, scenarios M1 (latent loading constant), M3 (loading ~ time),
# M5 (loading ~ time + dose). Output: tables/V4_S6_mmm_se_probe.csv, ..._summary.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
suppressMessages({ library(MMLB); library(dplyr) })
RD <- file.path(REPO, "R/v4/rundir")
for (f in c("generation_model1.R", "generation_model3.R", "generation_model5.R", "generate_data.R", "make_all_data.R")) source(file.path(RD, f))
FORM <- list(`1` = list(m = binY ~ time * trt, lv = ~1, t = ~1), `3` = list(m = binY ~ time * trt, lv = ~time, t = ~1),
             `5` = list(m = binY ~ time * trt + dose, lv = ~time + dose, t = ~trt + steroids))
cells <- data.table(model = c(1L, 3L, 5L), N = 145L, ntime = 10L, label = "manual")
one <- function(model, N, ntime, label, i) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 1000L * model + i)
  d <- tryCatch(as.data.frame(make_all_data(N, ntime, i, model, label)), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  F <- FORM[[as.character(model)]]
  fit <- tryCatch(mm(F$m, lv.formula = F$lv, t.formula = F$t, data = d, id = id, q = 20), error = function(e) NULL)
  if (is.null(fit)) return(data.table(model = model, i = i, ok = FALSE))
  k <- which(names(fit$beta) == "time:trt"); if (!length(k)) k <- grep("time:trt", rownames(fit$mod.cov))[1]
  data.table(model = model, i = i, ok = TRUE, est = unname(fit$beta[k]), se_mod = sqrt(fit$mod.cov[k, k]), se_rob = sqrt(fit$rob.cov[k, k]),
             se_summary = summary(fit)$mean.table[k, 2])
}
jobs <- cells[, .(i = 1:200), by = .(model, N, ntime, label)]
R <- rbindlist(parallel::mclapply(seq_len(nrow(jobs)), function(r) with(jobs[r], one(model, N, ntime, label, i)), mc.cores = 30, mc.preschedule = FALSE), fill = TRUE)
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[method == "mmm" & dgp == "GenBinaryY" & N == 145 & ntime == 10 & label == "manual" & term == "time:trt" & sim_id <= 200]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
R <- merge(R, V[, .(model, i = sim_id, est_stored = estimate, se_stored = se)], by = c("model", "i"), all.x = TRUE)
fwrite(R, file.path(TABDIR, "V4_S6_mmm_se_probe.csv"))
source(file.path(REPO, "R/s6/s6_mmm_se_probe_sum.R"))
