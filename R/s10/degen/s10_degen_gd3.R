# =============================================================================
# s10_degen_gd3.R -- S10 "Degeneri", gate GD3: on a stratified sample (seed 20261102; 2 eligible datasets per design
# point and arm, i.e. 4 per point, 9,504 in all) every dataset is regenerated and checked with the LP (lp_sep) AND with
# detectseparation 0.4.0, whatever the LP says; agreement must be total. Eligibility (generated, estimable design, outcome
# not constant) is read from the stored recs (S7d, S7e) so that the sample does not depend on the main run.
# Also, as a descriptive check that does not use linear programming: a plain glm fit (maxit 100) with its warnings.
# Outputs: results_v4/S10/degen/gd3/gd3_rows.rds, tables/V7_S10_degen_gd3.csv
# =============================================================================
source(path.expand("~/MMM/R/s10/degen/s10_degen_common.R")); setDTthreads(1L)
NC <- as.integer(Sys.getenv("GD3_CORES", "4"))
TK <- dg_tasks()
R <- rbindlist(parallel::mclapply(seq_len(nrow(TK)), function(r) { t <- TK[r]
  readRDS(t$rec_file)$recs[, .(src = t$src, dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label, chunk = t$chunk,
                               sim_id = as.integer(sim_id), ok = generated & !design_degenerate & sum_y > 0 & sum_y < n_obs)] }, mc.cores = NC))
R <- R[ok %in% TRUE]; setorder(R, src, dgp, scen, N, ntime, label, sim_id)
set.seed(20261102L)
SMP <- R[, .(sim_id = if (.N <= 2L) sim_id else sort(sample(sim_id, 2L))), by = .(src, dgp, scen, N, ntime, label)]
cat("sampled", nrow(SMP), "datasets;", uniqueN(SMP[, .(dgp, scen, N, ntime)]), "points\n")
SMP[, chunk := (sim_id - 1L) %/% 100L + 1L]
SMP <- merge(SMP, TK[, .(src, dgp, scen, N, ntime, label, chunk, model, rec_file, out)], by = c("src", "dgp", "scen", "N", "ntime", "label", "chunk"))
for (g in unique(SMP$dgp)) invisible(dg_env(g))
glm_chk <- function(t, d) { w <- character(0)
  f <- withCallingHandlers(tryCatch(glm(FORM_OF(t$scen), data = d, family = binomial, control = glm.control(maxit = 100)), error = function(e) NULL),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  if (is.null(f)) return(list(glm_ok = FALSE))
  cf <- coef(f); s <- summary(f)$coefficients; k <- "time:trt"
  list(glm_ok = TRUE, glm_conv = f$converged, glm_w01 = any(grepl("numerically 0 or 1", w)), glm_wnc = any(grepl("did not converge", w)),
       glm_maxabs = max(abs(cf), na.rm = TRUE), glm_b = unname(cf[k]), glm_se = if (k %in% rownames(s)) s[k, 2] else NA_real_) }
OUT <- rbindlist(parallel::mclapply(seq_len(nrow(SMP)), function(r) { t <- SMP[r]; gen <- dg_gen(t)
  x <- dg_one(gen, t, t$sim_id, ds_always = TRUE)
  set_rep_seed(t$model, t$sim_id); d <- as.data.frame(gen(t$N, t$ntime, t$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  c(list(src = t$src, dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label), x, glm_chk(t, d)) }, mc.cores = NC, mc.preschedule = FALSE), fill = TRUE)
saveRDS(OUT, file.path(OUT_DG, "gd3/gd3_rows.rds"))
pp <- OUT[, .N, by = .(dgp, scen, N, ntime)]
G3 <- data.table(check = c("datasets in the sample", "design points covered", "design points with fewer than 2 datasets", "eligible after regeneration",
                           "LP status not 0", "detectseparation errors", "LP separated", "detectseparation separated", "LP and detectseparation agree",
                           "time:trt infinite (detectseparation)", "glm: warning fitted 0/1 or not converged, among LP separated",
                           "glm: warning fitted 0/1 or not converged, among LP not separated", "glm: |b time:trt| > 10, among LP separated", "glm: |b time:trt| > 10, among LP not separated"),
                 value = c(nrow(OUT), nrow(pp), sum(pp$N < 2L), OUT[, sum(elig)], OUT[, sum(is.na(lp))], OUT[, sum(!is.na(ds_err))], OUT[, sum(lp %in% TRUE)],
                           OUT[, sum(ds_outcome %in% TRUE)], OUT[, sum(identical_ <- (lp %in% TRUE) == (ds_outcome %in% TRUE) & !is.na(lp) & !is.na(ds_outcome))],
                           OUT[, sum(grepl("(^|;)time:trt(;|$)", ds_inf_terms))], OUT[lp %in% TRUE, sum(glm_w01 | glm_wnc | !glm_conv, na.rm = TRUE)],
                           OUT[lp %in% FALSE, sum(glm_w01 | glm_wnc | !glm_conv, na.rm = TRUE)], OUT[lp %in% TRUE, sum(abs(glm_b) > 10, na.rm = TRUE)], OUT[lp %in% FALSE, sum(abs(glm_b) > 10, na.rm = TRUE)]))
G3[, pass := NA]; G3[check == "LP and detectseparation agree", pass := value == nrow(OUT)]
G3[check == "datasets in the sample", pass := value >= 5000]; G3[check == "design points with fewer than 2 datasets", pass := value == 0]
fwrite(G3, "/home/user/MMM/tables/V7_S10_degen_gd3.csv"); print(G3)
