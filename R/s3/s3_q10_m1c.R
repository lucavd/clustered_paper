# S3 Q10: calibration of scenario M1c (D-14), verified before the run.
source(path.expand("~/MMM/R/s3/s3_common.R"))
library(data.table)
E <- lapply(setNames(DGPS, DGPS), function(d) load_env(d, RUNDIR_V4)); P <- E[["GenBinaryY"]]
coh <- fread(file.path(REPO, "tables/V4_S2_lus_vs_dgm.csv"))
pre <- c(`2` = 0.2427, `3` = 0.1960, `4` = 0.1495, `6` = 0.0569)
# Q10.5a constants vs the cohort table
q105a <- data.table(check = c("PREV == prevalence_lus_severe", "VAR == var_id_modelA", "SIGMA_GBY == sqrt(VAR)", "GAMMA_GBY == 1",
                              "ZEGER_F == sqrt(1 + c^2 VAR)", "AR1_RHO == 0.6", "AR1_SIGMA2 == 0.5"),
  pass = c(P$M1C_PREV == coh[item == "prevalence_lus_severe", value], P$M1C_VAR == coh[item == "var_id_modelA", value],
           P$M1C_SIGMA_GBY == sqrt(P$M1C_VAR), P$M1C_GAMMA_GBY == 1,
           P$M1C_ZEGER_F == sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * P$M1C_VAR), E[["GLMM-AR1"]]$AR1_RHO == 0.6, E[["GLMM-AR1"]]$AR1_SIGMA2 == 0.5))
print(q105a); cat(sprintf("constants: sigma_GBY %.6f | Zeger f %.6f\n", P$M1C_SIGMA_GBY, P$M1C_ZEGER_F))
# Q10.1
gh <- fastGHQuad::gaussHermiteData(80); z <- gh$x * sqrt(2); wz <- gh$w / sqrt(pi)
prev_glmm <- function(T, inter) { t <- 0:(T - 1); w <- P$m1c_w(T); b <- P$M1C_ZEGER_F * c(P$m1c_beta0(T), 0.1, 0.2, inter)
  pm <- function(trt) vapply(t, function(tt) sum(wz * plogis(b[1] + b[2] * tt + b[3] * trt + b[4] * tt * trt + sqrt(P$M1C_VAR) * z)), 0)
  sum(w * (0.5 * pm(0) + 0.5 * pm(1))) / sum(w) }
q101 <- rbindlist(lapply(c(2, 3, 4, 6), function(T) { b0 <- P$m1c_beta0(T)
  data.table(ntime = T, beta0 = b0, prereg = pre[[as.character(T)]], resid = P$m1c_prev_marginal(b0, T, 0) - P$M1C_PREV,
             prev_H1_marginal = P$m1c_prev_marginal(b0, T, 0.2), prev_glmm_H0 = prev_glmm(T, 0), prev_glmm_H1 = prev_glmm(T, 0.2)) }))
q101[, pass_101 := abs(beta0 - prereg) < 5e-4 & abs(resid) < 1e-10]
q101[, pass_102 := abs(resid) < 1e-8 & abs(prev_glmm_H0 - 0.589) <= 0.001]
print(q101)
# Q10.5b: arguments actually handed to the three engines (N = 45, T = 6, i = 1)
cap <- function(dgp, engine, lab) { env <- E[[dgp]]; got <- NULL; real <- if (engine == "GenBinaryY") MMLB::GenBinaryY else get(engine, envir = env)
  assign(engine, function(...) { got <<- list(...); real(...) }, envir = env); on.exit(if (engine == "GenBinaryY") rm("GenBinaryY", envir = env) else assign(engine, real, envir = env))
  set_rep_seed(1L, 1); gen_data(env, "M1c", 1L, 45, 6, 1, lab); got }
b6 <- c(P$m1c_beta0(6), 0.1, 0.2, 0.2)
g1 <- cap("GenBinaryY", "GenBinaryY", "mmlb"); g2 <- cap("GLMM-native", "gen_glmm_data", "mmlb"); g3 <- cap("GLMM-AR1", "gen_ar1_data", "mmlb"); g3c <- cap("GLMM-AR1", "calibrate_beta_c_ar1", "mmlb")
q105b <- data.table(check = c("GBY beta", "GBY sigma", "GBY gamma", "GBY lv/t formula ~1", "GLMM beta_c = f x marginal", "GLMM D", "GLMM re ~1",
                              "AR1 target", "AR1 D (calib, gen)", "AR1 rho", "AR1 sigma2_ar"),
  pass = c(identical(g1$beta, b6), identical(g1$sigma, sqrt(P$M1C_VAR)), identical(g1$gamma, 1),
           deparse(g1$lv.formula) == "~1" && deparse(g1$t.formula) == "~1",
           isTRUE(all.equal(g2$beta_c / b6, rep(P$M1C_ZEGER_F, 4), tolerance = 1e-14)), identical(g2$D, P$M1C_VAR), deparse(g2$re.formula) == "~1",
           identical(g3c$target_beta_m, b6), identical(g3c$D, P$M1C_VAR) && identical(g3$D, P$M1C_VAR), identical(g3$rho, 0.6), identical(g3$sigma2_ar, 0.5)))
print(q105b); q105a <- rbind(q105a, q105b)
# Q10.3 / Q10.4: replicates 1-1000 at N = 45
J <- CJ(dgp = DGPS, ntime = c(2, 3, 4, 6), label = c("mmlb", "manual"))
emp <- rbindlist(lapply(seq_len(nrow(J)), function(k) {
  j <- J[k]; env <- E[[j$dgp]]
  r <- parallel::mclapply(1:1000, function(i) { deg <- design_degenerate(1L, 45, j$ntime, i); set_rep_seed(1L, i)
    d <- tryCatch(gen_data(env, "M1c", 1L, 45, j$ntime, i, j$label), error = function(e) NULL)
    if (is.null(d)) return(c(i, deg, NA, NA, NA))
    me <- attr(d, "calib_max_error"); c(i, deg, nrow(d), sum(d$binY), if (is.null(me)) NA else me) }, mc.cores = 45)
  m <- do.call(rbind, r); n <- m[, 3]; y <- m[, 4]; ok <- !is.na(n)
  p <- sum(y[ok]) / sum(n[ok]); u <- (y[ok] - p * n[ok]) / mean(n[ok])
  data.table(j, generated = sum(ok), degenerate = sum(m[, 2] == 1), prev = p, mcse = sd(u) / sqrt(sum(ok)),
             max_calib_err_nondeg = if (j$dgp == "GLMM-AR1") max(m[ok & m[, 2] == 0, 5]) else NA_real_)
}))
emp[, expected := ifelse(dgp == "GLMM-native",
  mapply(function(T, l) prev_glmm(T, if (l == "mmlb") 0.2 else 0), ntime, label),
  mapply(function(T, l) P$m1c_prev_marginal(P$m1c_beta0(T), T, if (l == "mmlb") 0.2 else 0), ntime, label))]
emp[, z := (prev - expected) / mcse][, pass_103 := abs(z) <= 3]
print(emp)
cat(sprintf("Q10.1 %d/4 | Q10.2 %d/4 | Q10.3 %d/24 (max |z| %.2f) | Q10.4 AR1 max calib error (non-degenerate) %.2e | Q10.5a %d/%d\n",
  sum(q101$pass_101), sum(q101$pass_102), sum(emp$pass_103), max(abs(emp$z)), max(emp$max_calib_err_nondeg, na.rm = TRUE), sum(q105a$pass), nrow(q105a)))
fwrite(q101, file.path(REPO, "tables/V4_S3_m1c_calibration.csv")); fwrite(emp, file.path(REPO, "tables/V4_S3_m1c_prevalence.csv"))
fwrite(q105a, file.path(REPO, "tables/V4_S3_m1c_constants.csv"))
