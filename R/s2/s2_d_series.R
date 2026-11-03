# =============================================================================
# s2_d_series.R -- S2(d), Q4: realised series lengths. Regenerates only the
# design draw nclust = ifelse(runif(N) < 0.7, T, sample(seq_len(T-1), N, TRUE)),
# exactly as written in every generator, under the production seed and RNGkind.
# Controprova Q4.4: compares with nclust/trt of full make_all_data() datasets.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
grid <- CJ(model = 1:5, N = N_LEVELS, ntime = T_LEVELS)
one <- function(model, N, ntime) {
  full <- 0; single <- 0; nobs <- 0
  for (i in 1:1000) {
    set_rep_seed(model, i)
    nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))
    full <- full + sum(nclust == ntime); single <- single + sum(nclust == 1); nobs <- nobs + sum(nclust)
  }
  n <- 1000 * N
  data.table(model, N, ntime, n_subj = n, p_full = full / n, p_single = single / n,
             mean_obs = nobs / n)
}
S <- rbindlist(parallel::mcMap(one, grid$model, grid$N, grid$ntime, mc.cores = 40))
S[, mcse_full := mcse_p(0.7, n_subj)]
S[, z_full := (p_full - 0.7) / mcse_full]
S[, pass_q42 := abs(p_full - 0.7) <= 0.01]
fwrite(S, file.path(TABDIR, "V4_S2_series_length.csv"))
tot <- S[, .(p_full = sum(p_full * n_subj) / sum(n_subj), n = sum(n_subj))]
cat(sprintf("Q4.1 pooled p_full = %.5f (MCSE %.5f) | pass = %s\n", tot$p_full,
            mcse_p(0.7, tot$n), abs(tot$p_full - 0.7) <= 0.01))
cat(sprintf("Q4.2 cells out of 0.70 +/- 0.01: %d of %d\n", sum(!S$pass_q42), nrow(S)))
print(S[pass_q42 == FALSE])
cat("Q4.2 by N: exceedances\n"); print(S[, .(n_fail = sum(!pass_q42), max_abs_z = max(abs(z_full))), by = N])
T2 <- S[ntime == 2]
cat(sprintf("Q4.3 T=2: p_single in [%.4f, %.4f]; all within 0.30 +/- 0.01: %s\n",
            min(T2$p_single), max(T2$p_single), all(abs(T2$p_single - 0.3) <= 0.01)))
S[, ratio := mean_obs / ntime]
cat(sprintf("Q4.3 mean_obs/T in [%.4f, %.4f]; all within 0.85 +/- 0.02: %s\n",
            min(S$ratio), max(S$ratio), all(abs(S$ratio - 0.85) <= 0.02)))
cat("mean obs per subject by T (pooled):\n"); print(S[, .(mean_obs = mean(mean_obs), p_single = mean(p_single)), by = ntime])

# ---- Q4.4 controprova --------------------------------------------------------
fam <- setNames(lapply(DGPS, load_family, rundir = RUNDIR), DGPS)
cells <- list(c(10, 2), c(45, 6))
chk <- list(); k <- 0L
for (m in 1:5) for (cl in cells) for (i in 1:5) {
  N <- cl[1]; Tt <- cl[2]
  set_rep_seed(m, i)
  nc_ref <- ifelse(runif(N) < 0.7, Tt, sample(seq_len(Tt - 1), N, replace = TRUE))
  per <- list()
  for (g in DGPS) {
    set_rep_seed(m, i)
    d <- tryCatch(fam[[g]]$make_all_data(N, Tt, i, m, "mmlb"), error = function(e) conditionMessage(e))
    if (is.character(d)) {
      k <- k + 1L
      chk[[k]] <- data.table(dgp = g, model = m, N = N, ntime = Tt, i = i, nclust_identical = NA, gen_error = d)
      per[[g]] <- NULL; next
    }
    d <- d[order(d$id, d$time), ]
    nc <- as.integer(table(factor(d$id, levels = seq_len(N))))
    tr <- tapply(d$trt, factor(d$id, levels = seq_len(N)), `[`, 1)
    per[[g]] <- list(nc = nc, tr = as.numeric(tr))
    k <- k + 1L
    chk[[k]] <- data.table(dgp = g, model = m, N = N, ntime = Tt, i = i,
                           nclust_identical = identical(as.integer(nc_ref), nc), gen_error = NA_character_)
  }
  k <- k + 1L
  ok <- names(per)
  same <- if (length(ok) >= 2) all(vapply(ok[-1], function(g) identical(per[[ok[1]]]$nc, per[[g]]$nc) &&
            isTRUE(all.equal(per[[ok[1]]]$tr, per[[g]]$tr, tolerance = 0)), logical(1))) else NA
  chk[[k]] <- data.table(dgp = "across", model = m, N = N, ntime = Tt, i = i, nclust_identical = same,
                         gen_error = if (length(ok) < 3) paste("compared:", paste(ok, collapse = ",")) else NA_character_)
}
C <- rbindlist(chk, fill = TRUE)
cat("Q4.4 generation errors:\n"); print(C[!is.na(gen_error)])
fwrite(C, file.path(TABDIR, "V4_S2_series_length_check.csv"))
cat(sprintf("Q4.4 regenerated vs make_all_data: %d/%d identical\n", sum(C[dgp != "across"]$nclust_identical, na.rm = TRUE), nrow(C[dgp != "across" & is.na(gen_error)])))
cat(sprintf("Q4.4 across mechanisms (nclust and trt): %d/%d identical\n", sum(C[dgp == "across"]$nclust_identical, na.rm = TRUE), nrow(C[dgp == "across"])))
