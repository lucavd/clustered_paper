# =============================================================================
# generation_models_ar1.R
# Scenario generators M1-M5 for the THIRD data-generating process:
# GLMM with AR(1) residual dependence (see MMM_sim_manual/gen_ar1_data.R).
#
# Design mirrors R_manual/ (the GLMM-native DGP) exactly -- same covariates,
# same random-effects structures, same follow-up truncation -- so that the only
# thing differing across the three DGPs is the dependence mechanism.
#
# Unlike the other two DGPs, the conditional coefficients are NOT hard-coded:
# adding the AR(1) residual variance changes the conditional->marginal
# attenuation, so a fixed beta_c would no longer imply the intended marginal
# effect. calibrate_beta_c_ar1() solves for the beta_c whose implied marginal
# vector equals the target, guaranteeing all three DGPs are compared at the
# SAME marginal effect size.
#
# H0 branch differs from H1 in the TESTED parameter (time:trt) only.
# =============================================================================

# AR(1) process parameters, shared across scenarios.
# rho = 0.6 gives clearly detectable serial dependence without dominating the
# cluster-level random effect; sigma2_ar = 0.5 keeps the total latent variance
# (D + sigma2_ar) comparable to the other DGPs.
AR1_RHO       <- 0.6
AR1_SIGMA2    <- 0.5

# Target MARGINAL coefficients, identical to the other two DGPs.
.ar1_target <- function(model, methods) {
  tgt <- switch(
    as.character(model),
    # (Intercept), time, trt, time:trt
    "1" = c(-2, 0.1, 0.2, 0.2),
    # (Intercept), time, trt, age, time:trt  -- interaction LAST (R term order)
    "2" = c(-2, 0.1, 0.2, 0.01, 0.2),
    "3" = c(-2, 0.1, 0.2, 0.2),
    "4" = c(-2, 0.2, 0.2, 0.05),
    # (Intercept), time, trt, dose, time:trt -- interaction LAST
    "5" = c(-2, 0.1, 0.2, 0.05, 0.2),
    stop("Unknown model")
  )
  if (methods == "manual") tgt[length(tgt)] <- 0   # H0: tested parameter only
  tgt
}

.ar1_design <- function(N, ntime, extra = c("none", "age", "full")) {
  extra <- match.arg(extra)
  # 70% of patients with full follow-up, 30% truncated
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))
  id   <- rep(seq(N), nclust)
  trt  <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  d <- data.frame(id, time, trt)
  if (extra %in% c("age", "full")) d$age <- scale(rep(rnorm(N, 50, 10), nclust))
  if (extra == "full") {
    d$sex      <- rep(rbinom(N, 1, 0.5), nclust)
    d$steroids <- rep(sample(1:3, N, replace = TRUE), nclust)
    # dose VARIES within subject, so a random slope on it is identified;
    # age and sex are drawn once per subject and would not be.
    d$dose     <- as.numeric(scale(rnorm(nrow(d))))
  }
  d[order(d$id, d$time), ]
}

.ar1_generate <- function(mean.formula, re.formula, D, data, model, methods, sim_id) {
  target <- .ar1_target(model, methods)
  beta_c <- calibrate_beta_c_ar1(
    mean.formula  = mean.formula,
    re.formula    = re.formula,
    target_beta_m = target,
    D             = D,
    sigma2_ar     = AR1_SIGMA2,
    data          = data
  )
  out <- gen_ar1_data(
    mean.formula = mean.formula,
    re.formula   = re.formula,
    beta_c       = as.numeric(beta_c),
    D            = D,
    rho          = AR1_RHO,
    sigma2_ar    = AR1_SIGMA2,
    time         = "time",
    id           = "id",
    data         = data,
    Yname        = "binY"
  )
  attr(out, "beta_m_target") <- target
  dplyr::mutate(out, sim_id = sim_id)
}


generate_data_model1 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  .ar1_generate(~ time * trt, ~ 1, D = 1.0,
                data = .ar1_design(N, ntime, "none"),
                model = 1, methods = methods, sim_id = sim_id)
}

generate_data_model2 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  .ar1_generate(~ time * trt + age, ~ 1, D = 1.0,
                data = .ar1_design(N, ntime, "age"),
                model = 2, methods = methods, sim_id = sim_id)
}

generate_data_model3 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  # random intercept + random slope on TIME (varies within subject; `age` does not)
  D <- matrix(c(1.0, 0.01,
                0.01, 0.001), nrow = 2, byrow = TRUE)
  .ar1_generate(~ time * trt, ~ 1 + time, D = D,
                data = .ar1_design(N, ntime, "age"),
                model = 3, methods = methods, sim_id = sim_id)
}

generate_data_model4 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  D <- matrix(c(1.0, 0.01,
                0.01, 0.001), nrow = 2, byrow = TRUE)
  .ar1_generate(~ time * trt, ~ 1 + time, D = D,
                data = .ar1_design(N, ntime, "age"),
                model = 4, methods = methods, sim_id = sim_id)
}

generate_data_model5 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  D <- matrix(c(1.0,  0.01,  0.01,
                0.01, 0.001, 0.0,
                0.01, 0.0,   0.01), nrow = 3, byrow = TRUE)
  .ar1_generate(~ time * trt + dose, ~ 1 + time + dose, D = D,
                data = .ar1_design(N, ntime, "full"),
                model = 5, methods = methods, sim_id = sim_id)
}
