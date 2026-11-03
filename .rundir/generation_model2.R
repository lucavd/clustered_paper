generate_data_model2 <- function(
  N,
  ntime,
  sim_id,
  methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # NB: R orders the terms of ~time*trt + age as
  #     (Intercept), time, trt, age, time:trt  -- the interaction comes LAST.
  # The previous vector c(-2, 0.1, 0.2, 0.2, 0.01) therefore set age = 0.2 and
  # time:trt = 0.01, while every downstream analysis scored M2 against
  # beta_m = 0.2. That transposition produced the spurious M2 anomalies
  # (coverage ~60%, MAPE ~118%, "power" == type I error). Fixed 2026-08-05.
  #
  # H0 branch differs from H1 in the tested parameter ONLY (time:trt = 0);
  # intercept and main effects are held at their H1 values so that type I error
  # and power are evaluated at the same outcome prevalence.
  beta_vec <- switch(
    methods,
    # (Intercept), time, trt, age, time:trt
    mmlb   = c(-2, 0.1, 0.2, 0.01, 0.2),
    manual = c(-2, 0.1, 0.2, 0.01, 0)
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

  beta_vec <- beta_vec

  sigma <- switch(
    methods,
    mmlb = 1,
    manual = 1)

  gamma <- switch(
    methods,
    mmlb = 1,
    manual = 1)

  GenBinaryY(
    mean.formula = ~ time*trt + age,
    lv.formula = ~1,
    t.formula = ~1,
    beta = beta_vec,

    sigma = sigma,
    gamma = gamma,
    id = id,
    data = data,
    q = 20,
    Yname = "binY"
  ) |>
    dplyr::mutate(
      sim_id = sim_id
    )
}
