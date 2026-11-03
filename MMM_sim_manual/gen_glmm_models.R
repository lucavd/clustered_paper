# =============================================================================
# gen_glmm_models.R
# Model-specific wrappers for GLMM-native data generation
#
# Analogous to R/generation_model1.R through R/generation_model5.R
# but using gen_glmm_data() instead of MMLB::GenBinaryY.
#
# Each function creates the covariate structure (with dropout) and calls
# gen_glmm_data() with model-specific parameters.
#
# methods = "alternative" → effect exists (for power/MAPE evaluation)
# methods = "null"        → no effect (for Type I error evaluation)
# =============================================================================

source("gen_glmm_data.R")


# =============================================================================
# Model 1: time*trt, random intercept
# MMLB analog: mean=~time*trt, lv=~1, t=~1, beta=(-2,0.1,0.2,0.2), sigma=1, gamma=1
# =============================================================================
gen_glmm_model1 <- function(N, ntime, sim_id,
                            methods = c("alternative", "null")) {
  methods <- match.arg(methods)
  set.seed(sim_id)

  # Conditional parameters (subject-specific)
  # Chosen so marginal ≈ MMLB's (-2, 0.1, 0.2, 0.2) after shrinkage
  beta_c <- switch(methods,
    alternative = c(-2.32, 0.116, 0.232, 0.232),
    null        = c(0, 0, 0, 0)
  )
  D <- 1.0  # random intercept variance

  # Dropout: 70% full, 30% truncated (same as current simulation)
  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id   <- rep(seq(N), nclust)
  trt  <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  data <- data.frame(id = id, time = time, trt = trt)
  data <- data[order(data$id, data$time), ]

  result <- gen_glmm_data(
    mean.formula = ~ time * trt,
    re.formula   = ~ 1,
    beta_c       = beta_c,
    D            = D,
    id           = id,
    data         = data,
    Yname        = "binY"
  )
  result$sim_id <- sim_id
  result
}


# =============================================================================
# Model 2: time*trt + age, random intercept
# MMLB analog: mean=~time*trt+age, lv=~1, t=~1, beta=(-2,0.1,0.2,0.2,0.01)
# =============================================================================
gen_glmm_model2 <- function(N, ntime, sim_id,
                            methods = c("alternative", "null")) {
  methods <- match.arg(methods)
  set.seed(sim_id)

  beta_c <- switch(methods,
    alternative = c(-2.32, 0.116, 0.232, 0.232, 0.012),
    null        = c(0, 0, 0, 0, 0)
  )
  D <- 1.0

  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id      <- rep(seq(N), nclust)
  trt     <- rep(rbinom(N, 1, 0.5), nclust)
  time    <- unlist(lapply(nclust, function(x) seq(x) - 1))
  age_raw <- rep(rnorm(N, mean = 50, sd = 10), nclust)
  age     <- as.numeric(scale(age_raw))
  data    <- data.frame(id = id, time = time, trt = trt, age = age)
  data    <- data[order(data$id, data$time), ]

  result <- gen_glmm_data(
    mean.formula = ~ time * trt + age,
    re.formula   = ~ 1,
    beta_c       = beta_c,
    D            = D,
    id           = id,
    data         = data,
    Yname        = "binY"
  )
  result$sim_id <- sim_id
  result
}


# =============================================================================
# Model 3: time*trt, random intercept + random slope on age
# MMLB analog: mean=~time*trt, lv=~age (sigma=c(1,0.01)), t=~1
# =============================================================================
gen_glmm_model3 <- function(N, ntime, sim_id,
                            methods = c("alternative", "null")) {
  methods <- match.arg(methods)
  set.seed(sim_id)

  beta_c <- switch(methods,
    alternative = c(-2.32, 0.116, 0.232, 0.232),
    null        = c(0, 0, 0, 0)
  )
  # Random intercept + random slope on age
  # Small random slope variance (mirroring MMLB's sigma = c(1, 0.01))
  D <- matrix(c(
    1.0,  0.01,
    0.01, 0.001
  ), nrow = 2, byrow = TRUE)

  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id      <- rep(seq(N), nclust)
  trt     <- rep(rbinom(N, 1, 0.5), nclust)
  time    <- unlist(lapply(nclust, function(x) seq(x) - 1))
  age_raw <- rep(rnorm(N, mean = 50, sd = 10), nclust)
  age     <- as.numeric(scale(age_raw))
  data    <- data.frame(id = id, time = time, trt = trt, age = age)
  data    <- data[order(data$id, data$time), ]

  result <- gen_glmm_data(
    mean.formula = ~ time * trt,
    re.formula   = ~ 1 + age,
    beta_c       = beta_c,
    D            = D,
    id           = id,
    data         = data,
    Yname        = "binY"
  )
  result$sim_id <- sim_id
  result
}


# =============================================================================
# Model 4: time*trt, random intercept + random slope on age
# MMLB analog: mean=~time*trt, lv=~age, t=~trt (adds trt-dependent transition)
# In GLMM DGP: same as Model 3 (no transition in GLMM).
# Different beta_c to match MMLB's (-2, 0.2, 0.2, 0.05).
# =============================================================================
gen_glmm_model4 <- function(N, ntime, sim_id,
                            methods = c("alternative", "null")) {
  methods <- match.arg(methods)
  set.seed(sim_id)

  beta_c <- switch(methods,
    alternative = c(-2.32, 0.232, 0.232, 0.058),
    null        = c(0, 0, 0, 0)
  )
  D <- matrix(c(
    1.0,  0.01,
    0.01, 0.001
  ), nrow = 2, byrow = TRUE)

  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id      <- rep(seq(N), nclust)
  trt     <- rep(rbinom(N, 1, 0.5), nclust)
  time    <- unlist(lapply(nclust, function(x) seq(x) - 1))
  age_raw <- rep(rnorm(N, mean = 50, sd = 10), nclust)
  age     <- as.numeric(scale(age_raw))
  data    <- data.frame(id = id, time = time, trt = trt, age = age)
  data    <- data[order(data$id, data$time), ]

  result <- gen_glmm_data(
    mean.formula = ~ time * trt,
    re.formula   = ~ 1 + age,
    beta_c       = beta_c,
    D            = D,
    id           = id,
    data         = data,
    Yname        = "binY"
  )
  result$sim_id <- sim_id
  result
}


# =============================================================================
# Model 5: time*trt, random intercept + random slopes on age and sex
# MMLB analog: mean=~time*trt, lv=~age+sex, t=~trt+steroids
# Most complex scenario.
# =============================================================================
gen_glmm_model5 <- function(N, ntime, sim_id,
                            methods = c("alternative", "null")) {
  methods <- match.arg(methods)
  set.seed(sim_id)

  beta_c <- switch(methods,
    alternative = c(-2.32, 0.116, 0.232, 0.232),
    null        = c(0, 0, 0, 0)
  )
  # Random intercept + random slopes on age and sex
  # Mirroring MMLB's sigma = c(1, 0.01, 0.1)
  D <- matrix(c(
    1.0,   0.01,  0.01,
    0.01,  0.001, 0.0,
    0.01,  0.0,   0.01
  ), nrow = 3, byrow = TRUE)

  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id       <- rep(seq(N), nclust)
  trt      <- rep(rbinom(N, 1, 0.5), nclust)
  time     <- unlist(lapply(nclust, function(x) seq(x) - 1))
  age_raw  <- rep(rnorm(N, mean = 50, sd = 10), nclust)
  age      <- as.numeric(scale(age_raw))
  sex      <- rep(rbinom(N, size = 1, prob = 0.5), nclust)
  steroids <- rep(sample(1:3, N, replace = TRUE), nclust)
  data     <- data.frame(id = id, time = time, trt = trt,
                         age = age, sex = sex, steroids = steroids)
  data     <- data[order(data$id, data$time), ]

  result <- gen_glmm_data(
    mean.formula = ~ time * trt,
    re.formula   = ~ 1 + age + sex,
    beta_c       = beta_c,
    D            = D,
    id           = id,
    data         = data,
    Yname        = "binY"
  )
  result$sim_id <- sim_id
  result
}


# =============================================================================
# Dispatcher (analogous to R/generate_data.R)
# =============================================================================
gen_glmm_dispatch <- function(model, N, ntime, sim_id,
                              methods = c("alternative", "null")) {
  methods <- match.arg(methods)
  switch(as.character(model),
    "1" = gen_glmm_model1(N, ntime, sim_id, methods),
    "2" = gen_glmm_model2(N, ntime, sim_id, methods),
    "3" = gen_glmm_model3(N, ntime, sim_id, methods),
    "4" = gen_glmm_model4(N, ntime, sim_id, methods),
    "5" = gen_glmm_model5(N, ntime, sim_id, methods),
    stop("Unknown model: ", model)
  )
}
