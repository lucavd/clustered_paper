# =============================================================================
# s2_c_truth.R -- S2(c), Q3: true values per mechanism x scenario x T (x arm).
# No model is fitted. GLMM-native marginal truth by exact Gauss-Hermite
# quadrature (protocol v1.2); GLMM-AR1 conditional truth from the generator's
# own calibration on regenerated designs.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
P <- fread(file.path(TABDIR, "V4_dgm_parameters.csv"))
getp <- function(g, m, lab, par) P[dgp == g & model == m & label == lab & parameter == par]$value
ARM <- c(mmlb = "H1", manual = "H0")

# ---- Q3.1 GenBinaryY -----------------------------------------------------------
GB <- rbindlist(lapply(1:5, function(m) rbindlist(lapply(names(ARM), function(lab) {
  b <- txt_num(getp("GenBinaryY", m, lab, "beta"))
  data.table(dgp = "GenBinaryY", model = m, label = lab, truth_marginal = tail(b, 1))
}))))
GB[, expected := fifelse(label == "manual", 0, fifelse(model == 4, 0.05, 0.2))]
cat(sprintf("Q3.1 GenBinaryY marginal truth equals 0.2/0.05/0: %d/%d\n", sum(GB$truth_marginal == GB$expected), nrow(GB)))

# ---- GLMM-native: marginal limit by quadrature ---------------------------------
gh_std <- function(q) { g <- fastGHQuad::gaussHermiteData(q); list(x = sqrt(2) * g$x, w = g$w / sqrt(pi)) }
# p(x) = E_b expit(x'beta_c + z'b), z'b ~ N(0, z'Dz)
marg_prob <- function(eta, v, gh) vapply(seq_along(eta), function(k) sum(gh$w * plogis(eta[k] + sqrt(v[k]) * gh$x)), 0)
# weighted logistic fit of the marginal mean model to exact probabilities (independence-GEE limit)
irls <- function(X, p, w) {
  beta <- rep(0, ncol(X))
  for (it in 1:500) {
    mu <- plogis(drop(X %*% beta)); W <- w * mu * (1 - mu)
    step <- solve(crossprod(X, X * W), crossprod(X, w * (p - mu)))
    beta <- beta + drop(step); if (max(abs(step)) < 1e-14) break
  }
  setNames(beta, colnames(X))
}
# design support: time x trt (x continuous covariate on GH nodes), with weights
design_grid <- function(m, T, weighted, gh_cov) {
  g <- CJ(time = 0:(T - 1), trt = 0:1)
  g[, w := if (weighted) fifelse(time == 0, 1, 0.7 + 0.3 * (T - 1 - time) / (T - 1)) * 0.5 else 1]
  if (m == 2 || m == 5) {
    cv <- data.table(cov = gh_cov$x, wc = gh_cov$w)
    g <- g[, as.list(cv), by = .(time, trt, w)][, w := w * wc][, wc := NULL]
  } else g[, cov := 0]
  g
}
native_truth <- function(m, lab, T, weighted, q = 80) {
  gh <- gh_std(q)
  bc <- txt_num(getp("GLMM-native", m, lab, "beta_c"))
  D  <- txt_mat(getp("GLMM-native", m, lab, "D"))
  g <- design_grid(m, T, weighted, gh)
  X <- switch(as.character(m),
    "2" = cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, age = g$cov, `time:trt` = g$time * g$trt),
    "5" = cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, dose = g$cov, `time:trt` = g$time * g$trt),
    cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, `time:trt` = g$time * g$trt))
  Z <- switch(as.character(m), "1" = matrix(1, nrow(X), 1), "2" = matrix(1, nrow(X), 1), "5" = cbind(1, g$time, g$cov), cbind(1, g$time))
  v <- rowSums((Z %*% D) * Z)
  p <- marg_prob(drop(X %*% bc), v, gh)
  b <- irls(X, p, g$w)
  unname(b["time:trt"])
}
NT <- rbindlist(lapply(1:5, function(m) rbindlist(lapply(names(ARM), function(lab) rbindlist(lapply(T_LEVELS, function(T) {
  bc <- txt_num(getp("GLMM-native", m, lab, "beta_c"))
  data.table(dgp = "GLMM-native", model = m, label = lab, ntime = T, truth_conditional = tail(bc, 1),
             truth_marginal_weighted = native_truth(m, lab, T, TRUE),
             truth_marginal_balanced = native_truth(m, lab, T, FALSE),
             quad_err_40v80 = abs(native_truth(m, lab, T, TRUE, 40) - native_truth(m, lab, T, TRUE, 80)))
}))))))
NT[, nominal_marginal := fifelse(label == "manual", 0, fifelse(model == 4, 0.05, 0.2))]
cat("Q3.2 GLMM-native conditional truth values:\n"); print(unique(NT[, .(model, label, truth_conditional)]))
cat(sprintf("Q3.4 max |40 - 80 nodes| = %.3g (pass <= 1e-10: %s)\n", max(NT$quad_err_40v80), max(NT$quad_err_40v80) <= 1e-10))

# Q3.3 reproduce F_truth_glmmnative.csv (balanced grid, D = 1 intercept only, M1 coefficients and M4 coefficients)
FT <- fread(file.path(TABDIR, "F_truth_glmmnative.csv"))
bal_D1 <- function(bc, T) {
  gh <- gh_std(80); g <- CJ(time = 0:(T - 1), trt = 0:1)
  X <- cbind(1, g$time, g$trt, g$time * g$trt); colnames(X) <- c("(Intercept)", "time", "trt", "time:trt")
  unname(irls(X, marg_prob(drop(X %*% bc), rep(1, nrow(X)), gh), rep(1, nrow(X)))["time:trt"])
}
FT[, mine_main := vapply(ntime, function(T) bal_D1(c(-2.32, 0.116, 0.232, 0.232), T), 0)]
FT[, mine_m4 := vapply(ntime, function(T) bal_D1(c(-2.32, 0.232, 0.232, 0.058), T), 0)]
FT[, `:=`(d_main = abs(mine_main - m_main), d_m4 = abs(mine_m4 - m_m4))]
print(FT[, .(ntime, m_main, mine_main, d_main, m_m4, mine_m4, d_m4)], digits = 12)
cat(sprintf("Q3.3 max abs diff vs F_truth_glmmnative = %.3g (pass <= 1e-10: %s)\n",
            max(FT$d_main, FT$d_m4), max(FT$d_main, FT$d_m4) <= 1e-10))

# ---- GLMM-AR1: conditional truth from the generator's calibration ---------------
fam <- load_family("GLMM-AR1", RUNDIR)
SPEC <- list("1" = list(f = ~ time * trt, re = ~ 1, ex = "none"),
             "2" = list(f = ~ time * trt + age, re = ~ 1, ex = "age"),
             "3" = list(f = ~ time * trt, re = ~ 1 + time, ex = "age"),
             "4" = list(f = ~ time * trt, re = ~ 1 + time, ex = "age"),
             "5" = list(f = ~ time * trt + dose, re = ~ 1 + time + dose, ex = "full"))
Dm <- function(m) txt_mat(getp("GLMM-AR1", m, "mmlb", "D"))
calib <- function(m, lab, N, T, i, rngkind = RNG_PROD) {
  RNGkind(rngkind[1], rngkind[2], rngkind[3]); set.seed(seed_of(m, i))
  sp <- SPEC[[as.character(m)]]
  d <- fam$.ar1_design(N, T, sp$ex)
  bc <- tryCatch(fam$calibrate_beta_c_ar1(mean.formula = sp$f, re.formula = sp$re,
                   target_beta_m = fam$.ar1_target(m, lab), D = Dm(m), sigma2_ar = fam$AR1_SIGMA2, data = d),
                 error = function(e) conditionMessage(e))
  if (is.character(bc)) return(list(val = NA_real_, err = bc, it = NA_integer_))
  list(val = tail(as.numeric(bc), 1), err = NA_character_, it = attr(bc, "iterations"))
}
grid <- CJ(model = 1:5, label = names(ARM), ntime = T_LEVELS, N = c(10, 25, 75, 145), i = 1:50)
res <- parallel::mcMap(function(m, lab, T, N, i) { r <- calib(m, lab, N, T, i); data.table(model = m, label = lab, ntime = T, N = N, i = i, beta_c = r$val, err = r$err, it = r$it) },
                       grid$model, grid$label, grid$ntime, grid$N, grid$i, mc.cores = 40)
AR <- rbindlist(res)
fwrite(AR, file.path(OUTDIR, "ar1_calibrated_by_replicate.csv"))
AS <- AR[, .(n_ok = sum(!is.na(beta_c)), n_err = sum(is.na(beta_c)), mean = mean(beta_c, na.rm = TRUE),
             sd = sd(beta_c, na.rm = TRUE), q025 = quantile(beta_c, 0.025, na.rm = TRUE), q975 = quantile(beta_c, 0.975, na.rm = TRUE),
             max_it = max(it, na.rm = TRUE)), by = .(model, label, ntime, N)]
AS[, target := fifelse(label == "manual", 0, fifelse(model == 4, 0.05, 0.2))]
fwrite(AS, file.path(TABDIR, "V4_S2_truth_ar1_conditional.csv"))
cat(sprintf("Q3.5a AR1 (model,label,T,N) with sd < 1e-8: %d of %d\n", sum(AS$sd < 1e-8, na.rm = TRUE), nrow(AS)))
cat("Q3.5a sd by T (max over model, label, N):\n"); print(AS[, .(max_sd = max(sd, na.rm = TRUE), min_sd = min(sd, na.rm = TRUE)), by = ntime])
cat("calibration errors:\n"); print(AR[!is.na(err), .N, by = .(model, ntime, N, err = substr(err, 1, 60))])
cat("AR1 conditional truth (mean over replicates), H1 and H0, N = 145:\n")
print(dcast(AS[N == 145], model + label ~ ntime, value.var = "mean"), digits = 4)

# Q3.5b reproduce F_truth_ar1.csv: its own procedure (Mersenne-Twister default, seeds rep 1..6, N 25/75/145).
FA <- fread(file.path(TABDIR, "F_truth_ar1.csv"))
rep_old <- FA[, {
  v <- vapply(1:6, function(r) calib(model, "mmlb", N, ntime, r, rngkind = c("Mersenne-Twister", "Inversion", "Rejection"))$val, 0)
  .(cond_old = cond_used, cond_mine = mean(v, na.rm = TRUE))
}, by = .(model, ntime, N)]
rep_old[, d := abs(cond_old - cond_mine)]
cat("Q3.5b |F_truth_ar1 - reproduction| by model (max):\n"); print(rep_old[, .(max_d = max(d), pass = max(d) <= 1e-8), by = model])
fwrite(rep_old, file.path(OUTDIR, "ar1_F_truth_reproduction.csv"))

# exploratory (not pre-registered): population independence-GEE limit implied by the AR1 mean beta_c
ar1_pop <- function(m, lab, T) {
  gh <- gh_std(80); bcv <- AR[model == m & label == lab & ntime == T & N == 145 & !is.na(beta_c)]
  if (!nrow(bcv)) return(NA_real_)
  RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(seed_of(m, 1))
  sp <- SPEC[[as.character(m)]]; d <- fam$.ar1_design(145, T, sp$ex)
  bc_full <- fam$calibrate_beta_c_ar1(mean.formula = sp$f, re.formula = sp$re, target_beta_m = fam$.ar1_target(m, lab),
                                      D = Dm(m), sigma2_ar = fam$AR1_SIGMA2, data = d)
  g <- design_grid(m, T, TRUE, gh)
  X <- switch(as.character(m),
    "2" = cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, age = g$cov, `time:trt` = g$time * g$trt),
    "5" = cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, dose = g$cov, `time:trt` = g$time * g$trt),
    cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, `time:trt` = g$time * g$trt))
  Z <- switch(as.character(m), "1" = matrix(1, nrow(X), 1), "2" = matrix(1, nrow(X), 1), "5" = cbind(1, g$time, g$cov), cbind(1, g$time))
  v <- rowSums((Z %*% Dm(m)) * Z) + fam$AR1_SIGMA2
  unname(irls(X, marg_prob(drop(X %*% as.numeric(bc_full)), v, gh), g$w)["time:trt"])
}
AP <- CJ(model = 1:5, label = names(ARM), ntime = T_LEVELS)
AP[, gee_limit_pop := mapply(ar1_pop, model, label, ntime)]
AP[, target := fifelse(label == "manual", 0, fifelse(model == 4, 0.05, 0.2))]
cat("exploratory: AR1 population GEE limit (beta_c of replicate 1 at N = 145) minus target:\n")
print(dcast(AP[, .(model, label, ntime, dev = gee_limit_pop - target)], model + label ~ ntime, value.var = "dev"), digits = 3)

# ---- assemble the v4 truth table ------------------------------------------------
TT <- rbindlist(list(
  GB[, .(dgp, model, label, ntime = NA_integer_, truth_marginal = truth_marginal, truth_marginal_balanced = truth_marginal,
         truth_conditional = NA_real_, conditional_note = "not a DGM parameter (GLMM target approximate)")],
  NT[, .(dgp, model, label, ntime, truth_marginal = truth_marginal_weighted, truth_marginal_balanced,
         truth_conditional, conditional_note = "fixed constant")],
  AS[N == 145, .(dgp = "GLMM-AR1", model, label, ntime, truth_marginal = target, truth_marginal_balanced = NA_real_,
                 truth_conditional = mean, conditional_note = sprintf("calibrated per replicate; mean of %d at N = 145, sd %.2g", n_ok, sd))]
), use.names = TRUE)
fwrite(TT, file.path(TABDIR, "V4_S2_truth_values.csv"))
fwrite(NT, file.path(TABDIR, "V4_S2_truth_glmmnative.csv"))
fwrite(AP, file.path(OUTDIR, "ar1_population_gee_limit_exploratory.csv"))
cat("GLMM-native marginal truth (weighted), by model x label x T:\n")
print(dcast(NT, model + label ~ ntime, value.var = "truth_marginal_weighted"), digits = 4)
cat("GLMM-native weighted minus balanced (max abs):", max(abs(NT$truth_marginal_weighted - NT$truth_marginal_balanced)), "\n")
