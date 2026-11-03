# =============================================================================
# gen_glmm_data.R
# Independent GLMM-native data generating process for multilevel binary data
#
# Generates data from a standard Generalized Linear Mixed Model:
#   logit(P(Y_ij = 1 | b_i)) = X_ij' beta_c + Z_ij' b_i
#   b_i ~ N(0, D)
#   Y_ij | b_i ~ indep Bernoulli(expit(X_ij' beta_c + Z_ij' b_i))
#
# NO dependency on MMLB package.
# NO transition/AR structure (conditional independence given random effects).
# Parameters are CONDITIONAL (subject-specific), not marginal.
#
# Also provides compute_marginal_truth() to derive the true marginal
# parameters that MMM/GEE target, via Gauss-Hermite quadrature.
# =============================================================================


#' Generate multilevel binary data from a standard GLMM
#'
#' @param mean.formula  One-sided formula for fixed effects (e.g., ~ time*trt)
#' @param re.formula    One-sided formula for random effects design matrix
#'                      (e.g., ~1 for random intercept, ~1+age for intercept+slope)
#' @param beta_c        Numeric vector of CONDITIONAL (subject-specific) coefficients.
#'                      Length must match ncol(model.matrix(mean.formula, data)).
#' @param D             Random effects covariance matrix. Can be:
#'                      - scalar (variance for random intercept)
#'                      - matrix (full covariance for multiple random effects)
#' @param id            Unquoted name of the cluster ID variable in data, OR a vector.
#' @param data          Data frame containing covariates and id variable.
#' @param Yname         Name for the generated binary outcome column (default "binY").
#'
#' @return data frame with all original columns plus:
#'   - Yname column: generated binary outcome
#'   - attr "beta_c": true conditional parameters
#'   - attr "beta_m_zeger": marginal parameters via Zeger approximation
#'   - attr "beta_m_exact": marginal parameters via GH quadrature (50 points)
#'   - attr "D": random effects covariance matrix
#'   - attr "dgp": "glmm" (to distinguish from GenBinaryY-generated data)
gen_glmm_data <- function(mean.formula, re.formula = ~1, beta_c, D,
                          id, data, Yname = "binY") {

  # --- Input validation ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.")

  id_vec <- if (is.character(id)) {
    data[[id]]
  } else {
    eval(substitute(id), data, parent.frame())
  }
  if (is.null(id_vec)) stop("Cannot find 'id' variable in data.")

  # --- Design matrices ---
  X <- model.matrix(mean.formula, model.frame(mean.formula, data))
  Z <- model.matrix(re.formula, model.frame(re.formula, data))

  if (ncol(X) != length(beta_c)) {
    stop(sprintf("beta_c length (%d) does not match design matrix columns (%d).",
                 length(beta_c), ncol(X)))
  }

  # --- Covariance matrix ---
  n_re <- ncol(Z)
  if (is.numeric(D) && length(D) == 1) {
    D <- matrix(D, 1, 1)
  }
  D <- as.matrix(D)
  if (nrow(D) != n_re || ncol(D) != n_re) {
    stop(sprintf("D must be %d x %d (matching re.formula columns), got %d x %d.",
                 n_re, n_re, nrow(D), ncol(D)))
  }

  # Verify D is positive semi-definite
  ev <- eigen(D, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < -1e-10)) {
    stop("D must be positive semi-definite.")
  }

  # --- Generate random effects ---
  unique_ids <- unique(id_vec)
  n_subj <- length(unique_ids)

  # Use MASS::mvrnorm for multivariate normal
  if (n_re == 1) {
    b_subj <- matrix(rnorm(n_subj, mean = 0, sd = sqrt(D[1, 1])), ncol = 1)
  } else {
    b_subj <- MASS::mvrnorm(n = n_subj, mu = rep(0, n_re), Sigma = D)
  }

  # Map subject-level random effects to observation level
  id_map <- match(id_vec, unique_ids)
  b_obs <- b_subj[id_map, , drop = FALSE]

  # --- Compute conditional probability ---
  # eta_ij = X_ij' beta_c + Z_ij' b_i
  eta <- as.numeric(X %*% beta_c) + rowSums(Z * b_obs)
  prob <- plogis(eta)

  # --- Generate binary outcome ---
  Y <- rbinom(length(prob), size = 1, prob = prob)

  # --- Compute true marginal parameters ---
  marginal <- compute_marginal_truth(X, beta_c, D, Z, q = 50)

  # --- Assemble output ---
  data[[Yname]] <- Y

  attr(data, "beta_c") <- beta_c
  attr(data, "beta_m_zeger") <- marginal$beta_m_zeger
  attr(data, "beta_m_exact") <- marginal$beta_m_exact
  attr(data, "D") <- D
  attr(data, "dgp") <- "glmm"
  attr(data, "shrinkage_factor") <- marginal$shrinkage_factor

  data
}


#' Compute true marginal parameters from conditional GLMM parameters
#'
#' For a logistic GLMM with random intercept:
#'   beta_m ≈ beta_c / sqrt(1 + c² D[1,1])         (Zeger approximation)
#'   beta_m_exact via Gauss-Hermite quadrature       (numerical)
#'
#' For random slopes: uses local linearization at covariate means.
#'
#' @param X       Fixed effects design matrix (n x p)
#' @param beta_c  Conditional coefficients (length p)
#' @param D       Random effects covariance matrix (q x q)
#' @param Z       Random effects design matrix (n x q)
#' @param q       Number of Gauss-Hermite quadrature points (default 50)
#'
#' @return List with:
#'   - beta_m_zeger: Zeger approximation of marginal parameters
#'   - beta_m_exact: GH-based marginal parameters
#'   - shrinkage_factor: Zeger shrinkage factor
compute_marginal_truth <- function(X, beta_c, D, Z = NULL, q = 50) {

  p <- length(beta_c)
  n_re <- ncol(D)

  # --- Zeger approximation ---
  # For random intercept: beta_m = beta_c / sqrt(1 + c^2 * D[1,1])
  # For random slopes: approximate using D[1,1] only (intercept variance dominates)
  c_const <- 16 * sqrt(3) / (15 * pi)
  shrink <- sqrt(1 + c_const^2 * D[1, 1])
  beta_m_zeger <- beta_c / shrink

  # --- Exact via Gauss-Hermite quadrature ---
  # For the random intercept model (n_re == 1):
  #   P(Y=1|X) = sum_k w_k * expit(X'beta_c + sqrt(D) * z_k)
  #   where (z_k, w_k) are GH nodes/weights for N(0,1)
  #
  # For random slopes (n_re > 1): use Monte Carlo integration
  # since multivariate GH is expensive

  if (n_re == 1) {
    beta_m_exact <- compute_marginal_exact_1d(X, beta_c, D[1, 1], q)
  } else {
    beta_m_exact <- compute_marginal_exact_mc(X, beta_c, D, Z, n_mc = 10000)
  }

  list(
    beta_m_zeger    = beta_m_zeger,
    beta_m_exact    = beta_m_exact,
    shrinkage_factor = shrink
  )
}


#' Exact marginal parameters via 1D Gauss-Hermite quadrature
#' (for random intercept models only)
compute_marginal_exact_1d <- function(X, beta_c, D_scalar, q = 50) {
  p <- length(beta_c)
  sigma <- sqrt(D_scalar)

  # Gauss-Hermite nodes and weights (for integral w.r.t. N(0,1))
  gh <- gauss_hermite_nodes(q)
  z_k <- gh$nodes   # scaled for N(0,1)
  w_k <- gh$weights  # scaled for N(0,1)

  # Strategy: compute marginal log-odds at reference points, extract coefficients
  # by finite differences.
  #
  # For a linear-in-X model, the marginal log-odds is approximately linear:
  #   logit(E[expit(X'beta_c + sigma*z)]) ≈ X' beta_m
  #
  # We compute beta_m by evaluating the marginal log-odds at p+1 carefully
  # chosen X values and solving the linear system.

  # Reference point: intercept = 1, all other covariates at zero
  x_ref <- rep(0, p)
  x_ref[1] <- 1  # intercept column is always 1 in model.matrix
  logit_p_ref <- compute_marginal_logodds(x_ref, beta_c, sigma, z_k, w_k)

  # For each coefficient j, perturb x_j by h and compute finite difference
  h <- 1e-4
  beta_m_exact <- numeric(p)
  beta_m_exact[1] <- logit_p_ref  # intercept = marginal log-odds at reference

  for (j in 2:p) {
    x_pert <- x_ref
    x_pert[j] <- x_pert[j] + h
    logit_p_pert <- compute_marginal_logodds(x_pert, beta_c, sigma, z_k, w_k)
    beta_m_exact[j] <- (logit_p_pert - logit_p_ref) / h
  }

  beta_m_exact
}


#' Exact marginal parameters via Monte Carlo integration
#' (for random slope models)
compute_marginal_exact_mc <- function(X, beta_c, D, Z, n_mc = 10000) {
  p <- length(beta_c)
  n_re <- ncol(D)

  # Draw MC samples of random effects
  if (n_re == 1) {
    b_mc <- matrix(rnorm(n_mc, 0, sqrt(D[1, 1])), ncol = 1)
  } else {
    b_mc <- MASS::mvrnorm(n_mc, mu = rep(0, n_re), Sigma = D)
  }

  # Same finite-difference approach but averaging over MC draws
  h <- 1e-4
  beta_m_exact <- numeric(p)

  # Use Z at the sample mean of covariates for the random effect part
  z_mean <- colMeans(Z)

  # Reference: intercept=1, all other covariates at zero
  x_ref <- rep(0, p)
  x_ref[1] <- 1  # intercept
  eta_ref <- sum(x_ref * beta_c)
  # For each MC draw, compute expit(eta_ref + z_mean' b)
  re_part <- as.numeric(b_mc %*% z_mean)
  p_ref <- mean(plogis(eta_ref + re_part))
  logit_ref <- qlogis(p_ref)
  beta_m_exact[1] <- logit_ref

  for (j in 2:p) {
    x_pert <- x_ref
    x_pert[j] <- h
    eta_pert <- sum(x_pert * beta_c)
    p_pert <- mean(plogis(eta_pert + re_part))
    logit_pert <- qlogis(p_pert)
    beta_m_exact[j] <- (logit_pert - logit_ref) / h
  }

  beta_m_exact
}


#' Compute marginal log-odds at a specific X value using GH quadrature
#'
#' logit(sum_k w_k * expit(x'beta_c + sigma * z_k))
compute_marginal_logodds <- function(x, beta_c, sigma, z_k, w_k) {
  eta_fixed <- sum(x * beta_c)
  # P(Y=1 | X=x) = sum_k w_k * expit(eta_fixed + sigma * z_k)
  p_marginal <- sum(w_k * plogis(eta_fixed + sigma * z_k))
  # Clamp to avoid Inf/-Inf
  p_marginal <- pmin(pmax(p_marginal, 1e-15), 1 - 1e-15)
  qlogis(p_marginal)
}


#' Gauss-Hermite quadrature nodes and weights (for N(0,1) integration)
#'
#' Returns nodes z_k and weights w_k such that:
#'   integral f(x) * phi(x) dx ≈ sum_k w_k * f(z_k)
#' where phi(x) is the standard normal density.
#'
#' No external package dependency — uses the Golub-Welsch algorithm.
gauss_hermite_nodes <- function(q) {
  # Golub-Welsch algorithm for Gauss-Hermite quadrature
  # Eigenvalues of the symmetric tridiagonal Jacobi matrix give nodes;
  # first component of eigenvectors gives weights.

  i <- seq_len(q - 1)
  # Recurrence coefficients for Hermite polynomials (physicists' convention)
  beta_vec <- sqrt(i / 2)

  # Symmetric tridiagonal matrix
  J <- matrix(0, q, q)
  diag(J) <- 0
  for (k in seq_len(q - 1)) {
    J[k, k + 1] <- beta_vec[k]
    J[k + 1, k] <- beta_vec[k]
  }

  eig <- eigen(J, symmetric = TRUE)
  nodes_raw <- eig$values
  weights_raw <- eig$vectors[1, ]^2 * sqrt(pi)

  # Sort by nodes
  ord <- order(nodes_raw)
  nodes_raw <- nodes_raw[ord]
  weights_raw <- weights_raw[ord]

  # Scale from physicists' Hermite (integral w.r.t. exp(-x^2))
  # to probabilists' (integral w.r.t. N(0,1) = phi(x)):
  # x_phys = x_prob * sqrt(2)
  # w_prob = w_phys / sqrt(pi) * sqrt(2) ... but we need to be careful.
  #
  # Standard: integral f(x) exp(-x^2) dx ≈ sum w_k f(x_k)
  # We want:  integral f(x) (1/sqrt(2pi)) exp(-x^2/2) dx ≈ sum w'_k f(x'_k)
  #
  # Substitution: x_phys = x_prob / sqrt(2)
  # dx_phys = dx_prob / sqrt(2)
  # exp(-x_phys^2) = exp(-x_prob^2/2)
  #
  # So: integral f(x_prob) * (1/sqrt(2pi)) * exp(-x_prob^2/2) dx_prob
  #   = integral f(sqrt(2)*x_phys) * (1/sqrt(2pi)) * exp(-x_phys^2) * sqrt(2) dx_phys
  #   ≈ sum_k w_k * f(sqrt(2)*x_phys_k) * sqrt(2) / sqrt(2pi)
  #   = sum_k (w_k / sqrt(pi)) * f(sqrt(2)*x_phys_k)

  nodes <- nodes_raw * sqrt(2)
  weights <- weights_raw / sqrt(pi)

  list(nodes = nodes, weights = weights)
}
