# =============================================================================
# m1c_params.R -- scenario M1c, anchored to the LUS cohort (D-13, D-14).
# Targets from tables/V4_S2_lus_vs_dgm.csv: realised marginal prevalence of
# lus_severe and latent patient-intercept variance (model A). The marginal
# intercept is calibrated per T on the H0 arm so that the EXPECTED realised
# prevalence over the truncated design equals the cohort prevalence; the same
# intercept is used under H1 (parity: arms differ in the tested coefficient only).
# Deterministic: no random numbers are consumed here.
# =============================================================================
M1C_PREV  <- 0.594512195121951     # prevalence_lus_severe
M1C_VAR   <- 3.54482462536461      # var_id_modelA
M1C_SIGMA_GBY <- sqrt(M1C_VAR)     # GenBinaryY latent loading (lv.formula = ~1)
M1C_GAMMA_GBY <- 1                 # GenBinaryY transition, as in M1
M1C_ZEGER_F   <- sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * M1C_VAR)  # GLMM-native, M1 convention

# expected number of observations at time t = 0..T-1 per subject:
# full follow-up with prob 0.7, otherwise uniform length on 1..T-1
m1c_w <- function(ntime) {
  t <- 0:(ntime - 1)
  0.7 + 0.3 * pmax(ntime - 1 - t, 0) / (ntime - 1)
}
m1c_prev_marginal <- function(b0, ntime, inter = 0) {
  t <- 0:(ntime - 1); w <- m1c_w(ntime)
  sum(w * (0.5 * plogis(b0 + 0.1 * t) + 0.5 * plogis(b0 + 0.1 * t + 0.2 + inter * t))) / sum(w)
}
m1c_beta0 <- function(ntime) {
  uniroot(function(b) m1c_prev_marginal(b, ntime, 0) - M1C_PREV,
          c(-5, 5), tol = 1e-14)$root
}
# (Intercept), time, trt, time:trt -- marginal
m1c_beta_marginal <- function(ntime, methods) {
  c(m1c_beta0(ntime), 0.1, 0.2, if (methods == "mmlb") 0.2 else 0)
}
