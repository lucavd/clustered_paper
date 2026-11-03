# =============================================================================
# gen_copula_data.R -- FOURTH data-generating family of the study (S7e, D-42):
# Gaussian copula with logistic margins ("dichotomised Gaussian"; Emrich and Piedmonte 1991;
# Song 2000; Song, Li and Yuan 2009; Masarotto and Varin 2012).
#
#   P(Y_ij = 1) = expit(X_ij' beta)                         exactly, for every observation
#   Y_ij = 1{Z_ij <= qnorm(expit(X_ij' beta))},  Z_i ~ N(0, R(t_i))
#   R(t)[j, k] = floor + (1 - floor) * rho^|t_j - t_k|      (j != k), R[j, j] = 1
#
# The marginal mean model is the study's mean model with the study's marginal coefficients, so
# the true marginal coefficient of time:trt is beta itself on every realised design: no
# calibration, no projection. The dependence is defined on the normal-score scale and does not
# involve the logit link: it is not a logistic mixed model (no subject-level logit-scale effect,
# conditional model not logistic), not MMLB's latent-variable + transition model, and on the
# Pearson scale the within-subject correlation is neither exchangeable nor AR(1) (it depends on
# the two margins). Random numbers: one rnorm(n_obs) call after the design draws, transformed
# subject by subject with the Cholesky factor of R for that subject's visit times.
# =============================================================================
copula_R <- function(tt, rho, floor) {
  R <- floor + (1 - floor) * rho^abs(outer(tt, tt, "-"))
  diag(R) <- 1
  R
}
gen_copula_data <- function(mean.formula, beta, rho, floor = 0, time = "time", id = "id",
                            data, Yname = "binY") {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (rho < 0 || rho >= 1) stop("'rho' must be in [0, 1).")
  if (floor < 0 || floor >= 1) stop("'floor' must be in [0, 1).")
  X <- model.matrix(mean.formula, model.frame(mean.formula, data))
  if (ncol(X) != length(beta)) stop(sprintf("beta length (%d) does not match design (%d).", length(beta), ncol(X)))
  mu <- plogis(as.numeric(X %*% beta))
  idv <- data[[id]]; tv <- data[[time]]
  e <- rnorm(length(idv))
  z <- numeric(length(idv))
  ch <- list()
  for (u in unique(idv)) {
    k <- which(idv == u); tt <- tv[k]; key <- paste(tt, collapse = ",")
    L <- ch[[key]]
    if (is.null(L)) { L <- t(chol(copula_R(tt, rho, floor))); ch[[key]] <- L }
    z[k] <- as.numeric(L %*% e[k])
  }
  data[[Yname]] <- as.integer(z <= qnorm(mu))
  attr(data, "beta_m_exact") <- stats::setNames(beta, colnames(X))
  attr(data, "copula") <- c(rho = rho, floor = floor)
  attr(data, "dgp") <- "gauss_copula"
  data
}
