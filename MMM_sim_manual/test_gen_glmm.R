# =============================================================================
# test_gen_glmm.R
# Comprehensive validation tests for gen_glmm_data()
#
# Tests:
#   1. Basic functionality (output structure, reproducibility)
#   2. No MMLB dependency
#   3. Gauss-Hermite quadrature accuracy
#   4. Marginal prevalence (empirical vs theoretical)
#   5. ICC (empirical vs theoretical)
#   6. Conditional parameter recovery (GLMM fit)
#   7. Marginal parameter recovery (GEE fit)
#   8. Zeger approximation accuracy (vs exact GH)
#
# Run: Rscript test_gen_glmm.R
# =============================================================================

cat("=============================================================\n")
cat("  VALIDATION TESTS FOR gen_glmm_data()\n")
cat("=============================================================\n\n")

source("gen_glmm_data.R")
source("gen_glmm_models.R")

n_pass <- 0
n_fail <- 0
n_test <- 0

check <- function(description, condition) {
  n_test <<- n_test + 1
  if (isTRUE(condition)) {
    n_pass <<- n_pass + 1
    cat(sprintf("  [PASS] %s\n", description))
  } else {
    n_fail <<- n_fail + 1
    cat(sprintf("  [FAIL] %s\n", description))
  }
}


# =============================================================================
# TEST 1: Basic functionality
# =============================================================================
cat("\n--- Test 1: Basic functionality ---\n")

d1 <- gen_glmm_model1(N = 50, ntime = 6, sim_id = 42, methods = "alternative")
check("Returns a data frame", is.data.frame(d1))
check("Has 'binY' column", "binY" %in% names(d1))
check("Has 'id' column", "id" %in% names(d1))
check("Has 'time' column", "time" %in% names(d1))
check("Has 'trt' column", "trt" %in% names(d1))
check("Has 'sim_id' column", "sim_id" %in% names(d1))
check("binY is binary (0/1)", all(d1$binY %in% c(0, 1)))
check("Has attr 'beta_c'", !is.null(attr(d1, "beta_c")))
check("Has attr 'beta_m_zeger'", !is.null(attr(d1, "beta_m_zeger")))
check("Has attr 'beta_m_exact'", !is.null(attr(d1, "beta_m_exact")))
check("Has attr 'D'", !is.null(attr(d1, "D")))
check("Has attr 'dgp' == 'glmm'", identical(attr(d1, "dgp"), "glmm"))
check("50 unique subjects", length(unique(d1$id)) == 50)


# =============================================================================
# TEST 2: Reproducibility
# =============================================================================
cat("\n--- Test 2: Reproducibility ---\n")

d2a <- gen_glmm_model1(N = 30, ntime = 4, sim_id = 123, methods = "alternative")
d2b <- gen_glmm_model1(N = 30, ntime = 4, sim_id = 123, methods = "alternative")
check("Same sim_id → same data", identical(d2a, d2b))

d2c <- gen_glmm_model1(N = 30, ntime = 4, sim_id = 124, methods = "alternative")
check("Different sim_id → different data", !identical(d2a$binY, d2c$binY))


# =============================================================================
# TEST 3: No MMLB dependency
# =============================================================================
cat("\n--- Test 3: No MMLB dependency ---\n")

# Check that gen_glmm_data.R does not load or call MMLB
src_text <- readLines("gen_glmm_data.R")
check("gen_glmm_data.R does not load MMLB",
      !any(grepl("library\\(MMLB\\)|require\\(MMLB\\)|MMLB::", src_text)))


# =============================================================================
# TEST 4: Gauss-Hermite quadrature accuracy
# =============================================================================
cat("\n--- Test 4: Gauss-Hermite quadrature accuracy ---\n")

# Test: integral of phi(x) dx from -Inf to Inf should equal 1
gh50 <- gauss_hermite_nodes(50)
int_1 <- sum(gh50$weights)
check(sprintf("GH integral of N(0,1) density = 1 (got %.6f)", int_1),
      abs(int_1 - 1) < 1e-6)

# Test: E[X] = 0
int_x <- sum(gh50$weights * gh50$nodes)
check(sprintf("GH E[X] for N(0,1) = 0 (got %.6f)", int_x),
      abs(int_x) < 1e-6)

# Test: E[X^2] = 1
int_x2 <- sum(gh50$weights * gh50$nodes^2)
check(sprintf("GH E[X^2] for N(0,1) = 1 (got %.6f)", int_x2),
      abs(int_x2 - 1) < 1e-6)

# Test: E[expit(1 + z)] should match numerical integration
gh_test <- sum(gh50$weights * plogis(1 + gh50$nodes))
# Compare with R's integrate()
ref_val <- integrate(function(z) plogis(1 + z) * dnorm(z), -10, 10)$value
check(sprintf("GH E[expit(1+z)] = %.6f vs integrate = %.6f", gh_test, ref_val),
      abs(gh_test - ref_val) < 1e-6)


# =============================================================================
# TEST 5: Marginal prevalence (empirical vs theoretical)
# =============================================================================
cat("\n--- Test 5: Marginal prevalence ---\n")

# Generate many datasets and check that empirical prevalence matches theory
# For intercept-only model: P(Y=1) = E[expit(beta_c0 + sigma*z)]
set.seed(999)
beta_c0 <- -1.0
D_test <- 1.0
sigma_test <- sqrt(D_test)
n_rep <- 2000
prev_vec <- numeric(n_rep)

for (r in seq_len(n_rep)) {
  N_test <- 100
  nclust_test <- rep(1, N_test)  # 1 obs per subject (simplest case)
  id_test <- seq(N_test)
  data_test <- data.frame(id = id_test)

  d_tmp <- gen_glmm_data(
    mean.formula = ~ 1,
    re.formula = ~ 1,
    beta_c = beta_c0,
    D = D_test,
    id = id_test,
    data = data_test,
    Yname = "Y"
  )
  prev_vec[r] <- mean(d_tmp$Y)
}

empirical_prev <- mean(prev_vec)
# Theoretical: E[expit(beta_c0 + sigma*z)] via GH
theoretical_prev <- sum(gh50$weights * plogis(beta_c0 + sigma_test * gh50$nodes))

check(sprintf("Marginal prevalence: empirical=%.3f, theoretical=%.3f (diff=%.4f)",
              empirical_prev, theoretical_prev, abs(empirical_prev - theoretical_prev)),
      abs(empirical_prev - theoretical_prev) < 0.02)


# =============================================================================
# TEST 6: ICC (empirical vs theoretical)
# =============================================================================
cat("\n--- Test 6: Intraclass correlation ---\n")

# For a logistic GLMM with random intercept variance D:
# ICC ≈ D / (D + pi^2/3)
D_icc <- 1.0
theoretical_icc <- D_icc / (D_icc + pi^2 / 3)

set.seed(888)
icc_vec <- numeric(200)
for (r in seq_len(200)) {
  N_icc <- 80
  nclust_icc <- rep(10, N_icc)
  id_icc <- rep(seq(N_icc), nclust_icc)
  data_icc <- data.frame(id = id_icc)

  d_icc <- gen_glmm_data(
    mean.formula = ~ 1,
    re.formula = ~ 1,
    beta_c = 0,  # 50% prevalence
    D = D_icc,
    id = id_icc,
    data = data_icc,
    Yname = "Y"
  )

  # Compute empirical ICC via ANOVA-type estimator
  grp_means <- tapply(d_icc$Y, d_icc$id, mean)
  grand_mean <- mean(d_icc$Y)
  msb <- 10 * var(grp_means)
  msw <- mean(tapply(d_icc$Y, d_icc$id, function(x) var(x) * (length(x)-1) / length(x)))
  if (msb > 0) {
    icc_vec[r] <- (msb - msw) / (msb + (10 - 1) * msw)
  } else {
    icc_vec[r] <- 0
  }
}

empirical_icc <- mean(icc_vec, na.rm = TRUE)
check(sprintf("ICC: empirical=%.3f, theoretical=%.3f (diff=%.3f)",
              empirical_icc, theoretical_icc, abs(empirical_icc - theoretical_icc)),
      abs(empirical_icc - theoretical_icc) < 0.05)


# =============================================================================
# TEST 7: Conditional parameter recovery (GLMM)
# =============================================================================
cat("\n--- Test 7: Conditional parameter recovery (GLMM) ---\n")

if (requireNamespace("lme4", quietly = TRUE)) {
  set.seed(777)
  beta_c_true <- c(-1.0, 0.3)
  D_true <- 1.0
  n_recov <- 500
  est_mat <- matrix(NA, n_recov, 2)

  for (r in seq_len(n_recov)) {
    N_r <- 100
    nclust_r <- rep(8, N_r)
    id_r <- rep(seq(N_r), nclust_r)
    x_r <- rnorm(length(id_r))
    data_r <- data.frame(id = id_r, x = x_r)

    d_r <- gen_glmm_data(
      mean.formula = ~ x,
      re.formula = ~ 1,
      beta_c = beta_c_true,
      D = D_true,
      id = id_r,
      data = data_r,
      Yname = "Y"
    )

    fit_r <- tryCatch(
      lme4::glmer(Y ~ x + (1 | id), data = d_r, family = binomial,
                  control = lme4::glmerControl(optimizer = "bobyqa")),
      error = function(e) NULL
    )
    if (!is.null(fit_r)) {
      est_mat[r, ] <- lme4::fixef(fit_r)
    }
  }

  mean_est <- colMeans(est_mat, na.rm = TRUE)
  convergence_rate <- mean(!is.na(est_mat[, 1]))

  check(sprintf("GLMM recovers intercept: true=%.2f, mean_est=%.3f (bias=%.3f)",
                beta_c_true[1], mean_est[1], mean_est[1] - beta_c_true[1]),
        abs(mean_est[1] - beta_c_true[1]) < 0.1)

  check(sprintf("GLMM recovers slope: true=%.2f, mean_est=%.3f (bias=%.3f)",
                beta_c_true[2], mean_est[2], mean_est[2] - beta_c_true[2]),
        abs(mean_est[2] - beta_c_true[2]) < 0.1)

  check(sprintf("GLMM convergence rate: %.1f%%", convergence_rate * 100),
        convergence_rate > 0.90)
} else {
  cat("  [SKIP] lme4 not available\n")
}


# =============================================================================
# TEST 8: Marginal parameter recovery (GEE)
# =============================================================================
cat("\n--- Test 8: Marginal parameter recovery (GEE) ---\n")

if (requireNamespace("geepack", quietly = TRUE)) {
  set.seed(666)
  beta_c_true2 <- c(-1.0, 0.3)
  D_true2 <- 1.0
  n_recov2 <- 500

  # Compute true marginal parameters
  gh <- gauss_hermite_nodes(50)
  sigma_gh <- sqrt(D_true2)
  # For intercept: marginal logodds at x=0
  p_at_0 <- sum(gh$weights * plogis(beta_c_true2[1] + sigma_gh * gh$nodes))
  beta_m_0_exact <- qlogis(p_at_0)
  # For slope: finite difference
  h_fd <- 1e-4
  p_at_h <- sum(gh$weights * plogis(beta_c_true2[1] + beta_c_true2[2]*h_fd + sigma_gh * gh$nodes))
  beta_m_1_exact <- (qlogis(p_at_h) - beta_m_0_exact) / h_fd
  beta_m_true <- c(beta_m_0_exact, beta_m_1_exact)

  est_gee <- matrix(NA, n_recov2, 2)
  for (r in seq_len(n_recov2)) {
    N_g <- 100
    nclust_g <- rep(8, N_g)
    id_g <- rep(seq(N_g), nclust_g)
    x_g <- rnorm(length(id_g))
    data_g <- data.frame(id = id_g, x = x_g)

    d_g <- gen_glmm_data(
      mean.formula = ~ x,
      re.formula = ~ 1,
      beta_c = beta_c_true2,
      D = D_true2,
      id = id_g,
      data = data_g,
      Yname = "Y"
    )
    d_g$id <- factor(d_g$id)

    fit_g <- tryCatch(
      geepack::geeglm(Y ~ x, id = id, data = d_g, family = binomial,
                       corstr = "exchangeable"),
      error = function(e) NULL
    )
    if (!is.null(fit_g)) {
      est_gee[r, ] <- coef(fit_g)
    }
  }

  mean_gee <- colMeans(est_gee, na.rm = TRUE)
  convergence_gee <- mean(!is.na(est_gee[, 1]))

  check(sprintf("GEE recovers marginal intercept: true=%.3f, mean_est=%.3f (bias=%.3f)",
                beta_m_true[1], mean_gee[1], mean_gee[1] - beta_m_true[1]),
        abs(mean_gee[1] - beta_m_true[1]) < 0.1)

  check(sprintf("GEE recovers marginal slope: true=%.3f, mean_est=%.3f (bias=%.3f)",
                beta_m_true[2], mean_gee[2], mean_gee[2] - beta_m_true[2]),
        abs(mean_gee[2] - beta_m_true[2]) < 0.1)
} else {
  cat("  [SKIP] geepack not available\n")
}


# =============================================================================
# TEST 9: Zeger approximation accuracy
# =============================================================================
cat("\n--- Test 9: Zeger approximation vs exact GH ---\n")

d9 <- gen_glmm_model1(N = 50, ntime = 6, sim_id = 42, methods = "alternative")
beta_m_z <- attr(d9, "beta_m_zeger")
beta_m_e <- attr(d9, "beta_m_exact")
beta_c_9 <- attr(d9, "beta_c")

cat(sprintf("  beta_c:       %s\n", paste(sprintf("%.4f", beta_c_9), collapse = ", ")))
cat(sprintf("  beta_m_zeger: %s\n", paste(sprintf("%.4f", beta_m_z), collapse = ", ")))
cat(sprintf("  beta_m_exact: %s\n", paste(sprintf("%.4f", beta_m_e), collapse = ", ")))

max_diff <- max(abs(beta_m_z - beta_m_e))
check(sprintf("Zeger approx close to exact GH (max diff = %.5f)", max_diff),
      max_diff < 0.05)


# =============================================================================
# TEST 10: All 5 models generate successfully
# =============================================================================
cat("\n--- Test 10: All models generate successfully ---\n")

for (m in 1:5) {
  for (meth in c("alternative", "null")) {
    d_test <- tryCatch(
      gen_glmm_dispatch(m, N = 30, ntime = 4, sim_id = 1, methods = meth),
      error = function(e) NULL
    )
    check(sprintf("Model %d, methods='%s' generates OK", m, meth),
          !is.null(d_test) && is.data.frame(d_test) && "binY" %in% names(d_test))
  }
}


# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=============================================================\n")
cat(sprintf("  RESULTS: %d passed, %d failed, %d total\n", n_pass, n_fail, n_test))
cat("=============================================================\n")

if (n_fail > 0) {
  cat("\n  WARNING: Some tests failed! Review output above.\n")
} else {
  cat("\n  All tests passed.\n")
}
