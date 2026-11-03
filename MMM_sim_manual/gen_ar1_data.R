# =============================================================================
# gen_ar1_data.R
# THIRD data-generating process: GLMM with AR(1) residual dependence.
#
# Rationale (referee comment, BMC MRM review 2026):
#   The study previously used two DGPs -- GenBinaryY (MMLB's own mechanism,
#   "home turf" for MMM) and gen_glmm_data (a standard GLMM, "home turf" for
#   GLMM). With only those two, the question "is MMM's advantage an artefact of
#   the design?" cannot be answered: whichever method wins, it wins on a DGP
#   built around one of the competitors. The referee also noted that in the
#   GLMM-native comparison the mechanism AND the target coefficient change
#   together, so estimand mismatch and model misspecification are confounded.
#
#   This DGP is deliberately NOBODY's home turf. It adds serial dependence that
#   a random intercept cannot absorb (violating the GLMM's conditional
#   independence assumption) but does so through a latent Gaussian AR(1) process
#   rather than MMLB's transition-model parameterisation (so it is not MMM's
#   mechanism either).
#
# Model:
#   logit(P(Y_ij = 1 | b_i, e_ij)) = X_ij' beta_c + Z_ij' b_i + e_ij
#   b_i  ~ N(0, D)                      cluster-level random effects
#   e_i  ~ N(0, sigma2_ar * R(rho))     within-cluster AR(1) residual process,
#                                       R(rho)[j,k] = rho^|t_j - t_k|
#   Y_ij | b_i, e_ij ~ indep Bernoulli
#
#   Setting rho = 0 and sigma2_ar = 0 recovers gen_glmm_data() exactly, which
#   makes this a strict generalisation and gives a continuous knob (rho) for
#   sensitivity analysis rather than a binary "misspecified / not" contrast.
#
# The true MARGINAL coefficients are obtained by integrating over BOTH random
# terms. Because b_i + e_ij is a sum of independent Gaussians, the marginal
# variance at observation level is D + sigma2_ar and the standard
# convolution result applies with that total variance -- see
# compute_marginal_truth_ar1() below, which does this exactly by quadrature
# rather than by the Zeger approximation.
#
# No dependency on MMLB.
# =============================================================================


#' Generate longitudinal binary data from a GLMM with AR(1) residual dependence
#'
#' @param mean.formula One-sided formula for fixed effects (e.g. ~ time*trt)
#' @param re.formula   One-sided formula for the random-effects design (e.g. ~1)
#' @param beta_c       CONDITIONAL (subject-specific) coefficients; length must
#'                     equal ncol(model.matrix(mean.formula, data))
#' @param D            Random-effects covariance: scalar or matrix
#' @param rho          AR(1) correlation of the latent residual process, |rho| < 1
#' @param sigma2_ar    Marginal variance of the latent residual process
#' @param time         Name of the time variable used for the AR(1) lag structure
#' @param id           Name of the cluster id variable, or a vector
#' @param data         Data frame of covariates
#' @param Yname        Name of the generated outcome column
#'
#' @return data with the outcome column added, carrying attributes
#'   beta_c, beta_m_exact, D, rho, sigma2_ar, dgp = "glmm_ar1"
gen_ar1_data <- function(mean.formula, re.formula = ~1, beta_c, D,
                         rho = 0.5, sigma2_ar = 0.5,
                         time = "time", id, data, Yname = "binY") {

  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (abs(rho) >= 1)        stop("'rho' must satisfy |rho| < 1.")
  if (sigma2_ar < 0)        stop("'sigma2_ar' must be non-negative.")

  id_vec <- if (is.character(id)) data[[id]] else eval(substitute(id), data, parent.frame())
  if (is.null(id_vec)) stop("Cannot find 'id' variable in data.")
  if (!time %in% names(data)) stop(sprintf("Time variable '%s' not found in data.", time))
  t_vec <- data[[time]]

  X <- model.matrix(mean.formula, model.frame(mean.formula, data))
  Z <- model.matrix(re.formula,   model.frame(re.formula,   data))

  if (ncol(X) != length(beta_c)) {
    stop(sprintf("beta_c length (%d) does not match design matrix columns (%d).",
                 length(beta_c), ncol(X)))
  }

  n_re <- ncol(Z)
  if (is.numeric(D) && length(D) == 1) D <- matrix(D, 1, 1)
  D <- as.matrix(D)
  if (nrow(D) != n_re || ncol(D) != n_re) {
    stop(sprintf("D must be %d x %d, got %d x %d.", n_re, n_re, nrow(D), ncol(D)))
  }
  if (any(eigen(D, symmetric = TRUE, only.values = TRUE)$values < -1e-10)) {
    stop("D must be positive semi-definite.")
  }

  # --- cluster-level random effects ---
  uid     <- unique(id_vec)
  n_subj  <- length(uid)
  b_subj  <- if (n_re == 1) matrix(rnorm(n_subj, 0, sqrt(D[1, 1])), ncol = 1)
             else MASS::mvrnorm(n_subj, rep(0, n_re), D)
  b_obs   <- b_subj[match(id_vec, uid), , drop = FALSE]

  # --- within-cluster AR(1) residual process, drawn per cluster on its own
  #     observed time points (handles unbalanced / truncated follow-up) ---
  e <- numeric(length(id_vec))
  if (sigma2_ar > 0) {
    for (u in uid) {
      idx <- which(id_vec == u)
      tt  <- t_vec[idx]
      R   <- rho ^ abs(outer(tt, tt, "-"))
      Sig <- sigma2_ar * R
      # jitter guards against numerical non-PSD from duplicated time points
      ch  <- tryCatch(chol(Sig + diag(1e-10, length(idx))),
                      error = function(err) chol(Sig + diag(1e-6, length(idx))))
      e[idx] <- as.numeric(crossprod(ch, rnorm(length(idx))))
    }
  }

  eta  <- as.numeric(X %*% beta_c) + rowSums(Z * b_obs) + e
  data[[Yname]] <- rbinom(length(eta), size = 1, prob = plogis(eta))

  beta_m <- compute_marginal_truth_ar1(X, beta_c, D, Z, sigma2_ar, q = 50)

  attr(data, "beta_c")       <- beta_c
  attr(data, "beta_m_exact") <- beta_m
  attr(data, "D")            <- D
  attr(data, "rho")          <- rho
  attr(data, "sigma2_ar")    <- sigma2_ar
  attr(data, "dgp")          <- "glmm_ar1"
  data
}


#' True marginal coefficients implied by conditional parameters under this DGP
#'
#' Integrates the conditional mean over the TOTAL observation-level latent
#' variance (random effects + AR(1) residual) by Gauss-Hermite quadrature, then
#' refits the marginal mean on the logit scale by weighted least squares over
#' the empirical covariate distribution. AR(1) correlation does not enter the
#' marginal MEAN -- only the dependence structure -- so rho is not an argument.
#'
#' @return numeric vector of marginal coefficients, named as the columns of X
compute_marginal_truth_ar1 <- function(X, beta_c, D, Z, sigma2_ar, q = 50) {
  gh <- fastGHQuad::gaussHermiteData(q)
  nodes <- gh$x * sqrt(2)
  wts   <- gh$w / sqrt(pi)

  # observation-level latent SD: z' D z (random effects) + sigma2_ar (AR process)
  v_re  <- rowSums((Z %*% D) * Z)
  sd_i  <- sqrt(v_re + sigma2_ar)

  eta_c <- as.numeric(X %*% beta_c)
  # marginal mean per observation: E_u[ expit(eta_c + sd_i * u) ]
  mu_m  <- vapply(seq_along(eta_c), function(i)
                    sum(wts * plogis(eta_c[i] + sd_i[i] * nodes)), numeric(1))

  # project onto the logit scale over the observed design
  eta_m <- qlogis(pmin(pmax(mu_m, 1e-12), 1 - 1e-12))
  fit   <- qr.solve(X, eta_m)
  stats::setNames(as.numeric(fit), colnames(X))
}


#' Calibrate beta_c so that the implied marginal coefficients match a target
#'
#' The three DGPs must be compared at the SAME marginal effect size, otherwise
#' differences in power/MAPE reflect effect size rather than the mechanism.
#' This solves for the conditional vector whose implied marginal vector equals
#' `target_beta_m`, by fixed-point iteration on the attenuation map.
calibrate_beta_c_ar1 <- function(mean.formula, re.formula = ~1, target_beta_m,
                                 D, sigma2_ar, data, tol = 1e-8, maxit = 200) {
  X <- model.matrix(mean.formula, model.frame(mean.formula, data))
  Z <- model.matrix(re.formula,   model.frame(re.formula,   data))
  if (is.numeric(D) && length(D) == 1) D <- matrix(D, 1, 1)

  beta_c <- target_beta_m               # start from the marginal values
  for (it in seq_len(maxit)) {
    implied <- compute_marginal_truth_ar1(X, beta_c, D, Z, sigma2_ar)
    err     <- target_beta_m - implied
    if (max(abs(err)) < tol) break
    beta_c <- beta_c + err              # damped-free update; contraction here
  }
  attr(beta_c, "iterations") <- it
  attr(beta_c, "max_error")  <- max(abs(err))
  beta_c
}
