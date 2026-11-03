# s6_mmm_se_probe_sum.R -- summary of tables/V4_S6_mmm_se_probe.csv (written by s6_mmm_se_probe.R).
# A replicate is usable, as in the study, when the estimate and the model-based SE are finite; the
# robust-SE columns use the same replicates (plus a finite robust SE).
source(path.expand("~/MMM/R/s6/s6_common.R"))
R <- fread(file.path(TABDIR, "V4_S6_mmm_se_probe.csv"))
z <- qnorm(0.975)
U <- R[ok == TRUE & is.finite(est) & is.finite(se_mod)]
S <- U[, .(n_fit = .N, n_rob = sum(is.finite(se_rob)),
           same_est = mean(abs(est - est_stored) < 1e-8), same_se = mean(abs(se_summary - se_stored) < 1e-8),
           emp_sd = sd(est), rms_se_mod = sqrt(mean(se_mod^2)), med_se_mod = median(se_mod),
           rms_se_rob = sqrt(mean(se_rob^2, na.rm = TRUE)), med_se_rob = median(se_rob, na.rm = TRUE),
           rej_mod = 100 * mean(abs(est / se_mod) > z), rej_rob = 100 * mean(abs(est / se_rob) > z, na.rm = TRUE)), by = model]
S <- merge(R[, .(n_rep = .N), by = model], S, by = "model")
S[, `:=`(ratio_mod = rms_se_mod / emp_sd, ratio_med_rob = med_se_rob / emp_sd, mcse_rej_mod = 100 * sqrt(rej_mod / 100 * (1 - rej_mod / 100) / n_fit),
         mcse_rej_rob = 100 * sqrt(rej_rob / 100 * (1 - rej_rob / 100) / n_rob))]
fwrite(S, file.path(TABDIR, "V4_S6_mmm_se_probe_summary.csv")); print(S, digits = 3)
