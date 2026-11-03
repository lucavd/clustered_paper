# s7e_gates_pre.R -- S7e gates checked BEFORE the run (docs/steps/S7e.md): GF2 (H0/H1 parity and coupling, same design as the
# other mechanisms) and GF3 (calibration of the generator against closed-form probabilities). No method is fitted to copula data.
source(path.expand("~/MMM/R/s7e/s7e_common.R"))
envC <- load_env_cop(); envN <- load_env("GLMM-native", RUNDIR_V4); envA <- load_env("GLMM-AR1", RUNDIR_V4)
FORM_MEAN <- function(m) switch(as.character(m), "2" = ~ time * trt + age, "5" = ~ time * trt + dose, ~ time * trt)
mu_of <- function(d, m, label) plogis(as.numeric(model.matrix(FORM_MEAN(m), d) %*% envC$.cop_target(m, label)))
gen_at <- function(env, m, N, Tt, i, label) { set_rep_seed(m, i); as.data.frame(env$make_all_data(N, Tt, i, m, label)) }
# ---- GF2 ----------------------------------------------------------------------------------------------------
G2 <- rbindlist(lapply(1:5, function(m) rbindlist(lapply(list(c(10L, 2L), c(25L, 6L)), function(nt) rbindlist(lapply(1:100, function(i) {
  d1 <- gen_at(envC, m, nt[1], nt[2], i, "mmlb"); d0 <- gen_at(envC, m, nt[1], nt[2], i, "manual")
  cols <- intersect(c("id", "time", "trt", "age", "dose"), names(d1))
  same_design <- identical(lapply(d1[cols], as.vector), lapply(d0[cols], as.vector))
  mu1 <- mu_of(d1, m, "mmlb"); mu0 <- mu_of(d0, m, "manual")
  mono <- all((mu1 > mu0 & d1$binY >= d0$binY) | (mu1 < mu0 & d1$binY <= d0$binY) | (mu1 == mu0 & d1$binY == d0$binY))
  # the GLMM-AR1 calibration fails on some non-estimable designs (generation error, D-11): compared only where it generates
  dn <- tryCatch(gen_at(envN, m, nt[1], nt[2], i, "manual"), error = function(e) NULL); da <- tryCatch(gen_at(envA, m, nt[1], nt[2], i, "manual"), error = function(e) NULL)
  same_N <- if (is.null(dn)) NA else identical(lapply(d0[cols], as.vector), lapply(dn[cols], as.vector))
  same_A <- if (is.null(da)) NA else identical(lapply(d0[cols], as.vector), lapply(da[cols], as.vector))
  data.table(model = m, N = nt[1], ntime = nt[2], sim_id = i, deg = design_degenerate(m, nt[1], nt[2], i), same_design, mono, same_N, same_A, n_diff_y = sum(d1$binY != d0$binY))
}))))))
g2 <- G2[, .(reps = .N, same_design = sum(same_design), mono = sum(mono), same_as_glmm_native = sum(same_N, na.rm = TRUE), same_as_glmm_ar1 = sum(same_A, na.rm = TRUE),
  gen_err_other = sum(is.na(same_N) | is.na(same_A)), gen_err_other_nondeg = sum((is.na(same_N) | is.na(same_A)) & !deg), mean_y_changed = mean(n_diff_y)), by = .(model, N, ntime)]
print(g2); fwrite(g2, file.path(TABDIR, "V4_S7e_GF2.csv"))
GF2_pass <- G2[, all(same_design & mono) && all(same_N[!is.na(same_N)]) && all(same_A[!is.na(same_A)]) && !any((is.na(same_N) | is.na(same_A)) & !deg)]
# ---- GF3 ----------------------------------------------------------------------------------------------------
pbv <- function(a, b, r) { if (!is.finite(a) || !is.finite(b)) return(NA_real_)   # P(Z1 <= a, Z2 <= b), corr r, by numerical integration
  integrate(function(x) dnorm(x) * pnorm((b - r * x) / sqrt(1 - r^2)), -Inf, a, rel.tol = 1e-10)$value }
NB <- 100000L; TB <- 10L; G3 <- list()
for (m in 1:5) {
  RNGkind(RNG_PROD[1], RNG_PROD[2], RNG_PROD[3]); set.seed(20261023L)
  d <- as.data.table(envC$make_all_data(NB, TB, 1L, m, "mmlb")); d[, mu := mu_of(as.data.frame(d), m, "mmlb")]
  if (m %in% c(1L, 3L, 4L)) {
    c1 <- d[, .(n = .N, obs = sum(binY), mu = mu[1]), by = .(time, trt)][, z := (obs - n * mu) / sqrt(n * mu * (1 - mu))]
    p <- merge(d[, .(id, t1 = time, trt, y1 = binY, mu1 = mu)], d[, .(id, t2 = time, y2 = binY, mu2 = mu)], by = "id", allow.cartesian = TRUE)[t2 > t1]
    c2 <- p[, .(n = .N, obs = sum(y1 * y2), mu1 = mu1[1], mu2 = mu2[1]), by = .(t1, t2, trt)]
    c2[, p11 := mapply(function(a, b, r) pbv(a, b, r), qnorm(mu1), qnorm(mu2), envC$COPULA_RHO^(t2 - t1))]
    c2[, z := (obs - n * p11) / sqrt(n * p11 * (1 - p11))]
    G3[[m]] <- data.table(model = m, check = c("marginal cells", "pair cells"), n_cells = c(nrow(c1), nrow(c2)),
                          max_abs_z = c(max(abs(c1$z)), max(abs(c2$z))), pass = c(all(abs(c1$z) < 4), all(abs(c2$z) < 4)))
  } else {
    X <- model.matrix(FORM_MEAN(m), as.data.frame(d)); y <- d$binY; bt <- envC$.cop_target(m, "mmlb")
    f <- glm.fit(X, y, family = binomial()); b <- f$coefficients; mu_h <- f$fitted.values
    A <- crossprod(X, X * (mu_h * (1 - mu_h))); U <- rowsum(X * (y - mu_h), d$id); Ai <- solve(A); V <- Ai %*% crossprod(U) %*% Ai
    z <- (b - bt) / sqrt(diag(V))
    G3[[m]] <- data.table(model = m, check = paste("GEE independence coef", colnames(X)), n_cells = 1L, max_abs_z = abs(z), pass = abs(z) < 4)
  }
}
G3 <- rbindlist(G3); print(G3, digits = 3); fwrite(G3, file.path(TABDIR, "V4_S7e_GF3.csv"))
cat(sprintf("GF2 %s | GF3 %s\n", if (GF2_pass) "PASS" else "FAIL", if (all(G3$pass)) "PASS" else "FAIL"))
