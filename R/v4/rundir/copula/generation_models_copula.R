# =============================================================================
# generation_models_copula.R -- scenarios M1-M5 for the Gaussian-copula family (S7e).
# Design draws are those of the GLMM families (R/v4/rundir/ar1/generation_models_ar1.R::.ar1_design),
# so under the study seed 20260805 + 1000 * model + i (L'Ecuyer-CMRG) each replicate has exactly the
# design (nclust, trt, age, dose, ...) of the other mechanisms; only the outcome draws differ.
# Marginal coefficients: the study's targets (H0 differs from H1 in time:trt only).
# The copula correlation is the same in the five scenarios; the scenarios differ in the mean model.
# COPULA_RHO / COPULA_FLOOR are set in copula_params.R (pre-registered in docs/steps/S7e.md).
# =============================================================================
.cop_target <- function(model, methods) {
  tgt <- switch(as.character(model),
    "1" = c(-2, 0.1, 0.2, 0.2), "2" = c(-2, 0.1, 0.2, 0.01, 0.2), "3" = c(-2, 0.1, 0.2, 0.2),
    "4" = c(-2, 0.2, 0.2, 0.05), "5" = c(-2, 0.1, 0.2, 0.05, 0.2), stop("Unknown model"))
  if (methods == "manual") tgt[length(tgt)] <- 0
  tgt
}
.cop_design <- function(N, ntime, extra = c("none", "age", "full")) {
  extra <- match.arg(extra)
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))
  id   <- rep(seq(N), nclust)
  trt  <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  d <- data.frame(id, time, trt)
  if (extra %in% c("age", "full")) d$age <- scale(rep(rnorm(N, 50, 10), nclust))
  if (extra == "full") {
    d$sex      <- rep(rbinom(N, 1, 0.5), nclust)
    d$steroids <- rep(sample(1:3, N, replace = TRUE), nclust)
    d$dose     <- as.numeric(scale(rnorm(nrow(d))))
  }
  d[order(d$id, d$time), ]
}
.cop_generate <- function(mean.formula, data, model, methods, sim_id) {
  out <- gen_copula_data(mean.formula, beta = .cop_target(model, methods), rho = COPULA_RHO,
                         floor = COPULA_FLOOR, time = "time", id = "id", data = data, Yname = "binY")
  dplyr::mutate(out, sim_id = sim_id)
}
generate_data_model1 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods); .cop_generate(~ time * trt, .cop_design(N, ntime, "none"), 1, methods, sim_id) }
generate_data_model2 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods); .cop_generate(~ time * trt + age, .cop_design(N, ntime, "age"), 2, methods, sim_id) }
generate_data_model3 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods); .cop_generate(~ time * trt, .cop_design(N, ntime, "age"), 3, methods, sim_id) }
generate_data_model4 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods); .cop_generate(~ time * trt, .cop_design(N, ntime, "age"), 4, methods, sim_id) }
generate_data_model5 <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods); .cop_generate(~ time * trt + dose, .cop_design(N, ntime, "full"), 5, methods, sim_id) }
