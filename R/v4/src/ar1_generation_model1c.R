# V4 (S3, D-13/D-14): scenario M1c for GLMM-AR1. Design from .ar1_design(, "none")
# as in M1; D = cohort variance with the AR(1) component on top (rho, sigma2_ar of
# M1); beta_c calibrated per replicate to the M1c marginal vector, as .ar1_generate.
generate_data_model1c <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  data   <- .ar1_design(N, ntime, "none")
  target <- m1c_beta_marginal(ntime, methods)
  beta_c <- calibrate_beta_c_ar1(
    mean.formula  = ~ time * trt,
    re.formula    = ~ 1,
    target_beta_m = target,
    D             = M1C_VAR,
    sigma2_ar     = AR1_SIGMA2,
    data          = data
  )
  out <- gen_ar1_data(
    mean.formula = ~ time * trt,
    re.formula   = ~ 1,
    beta_c       = as.numeric(beta_c),
    D            = M1C_VAR,
    rho          = AR1_RHO,
    sigma2_ar    = AR1_SIGMA2,
    time         = "time",
    id           = "id",
    data         = data,
    Yname        = "binY"
  )
  attr(out, "beta_m_target") <- target
  attr(out, "calib_max_error") <- attr(beta_c, "max_error")
  dplyr::mutate(out, sim_id = sim_id)
}
