generate_data_model2 <- function(
  N,
  ntime,
  sim_id,
  methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # Same transposition as R/generation_model2.R: R orders ~time*trt + age as
  # (Intercept), time, trt, age, time:trt -- the interaction LAST. The previous
  # vector therefore set age = 0.232 and time:trt = 0.012 while the analysis
  # scored M2 against 0.2. Fixed 2026-08-05, together with the H0 branch.
  beta_c <- switch(
    methods,
    # (Intercept), time, trt, age, time:trt -- conditional; marginal ≈ (-2, 0.1, 0.2, 0.01, 0.2)
    mmlb   = c(-2.32, 0.116, 0.232, 0.012, 0.232),
    manual = c(-2.32, 0.116, 0.232, 0.012, 0)
  )
  # 70% dei pazienti con ntime, 30% con un valore inferiore
  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id <- rep(seq(N), nclust)
  trt <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  age_raw <- rep(rnorm(N, mean = 50, sd = 10), nclust)
  age <- scale(age_raw)
  data <- data.frame(id, time, trt, age) |> arrange(id, time)

  D <- 1.0  # random intercept variance

  gen_glmm_data(
    mean.formula = ~ time*trt + age,
    re.formula = ~1,
    beta_c = beta_c,
    D = D,
    id = id,
    data = data,
    Yname = "binY"
  ) |>
    dplyr::mutate(
      sim_id = sim_id
    )
}
