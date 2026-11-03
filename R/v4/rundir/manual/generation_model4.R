generate_data_model4 <- function(
    N,
    ntime,
    sim_id,
    methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # M4 deliberately uses a small effect. H0 differs in the tested parameter only.
  beta_c <- switch(
    methods,
    # (Intercept), time, trt, time:trt -- conditional; marginal ≈ (-2, 0.2, 0.2, 0.05)
    mmlb   = c(-2.32, 0.232, 0.232, 0.058),
    manual = c(-2.32, 0.232, 0.232, 0)
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

  # Random intercept + random slope on TIME, which varies within subject.
  # An earlier version used `age`, drawn once per subject and repeated across
  # occasions: its random-effects design matrix had two proportional columns,
  # so the random slope WAS the random intercept rescaled -- not identified for
  # any number of occasions. Same reason `sex` was dropped from scenario 5.
  D <- matrix(c(
    1.0,  0.01,
    0.01, 0.001
  ), nrow = 2, byrow = TRUE)

  gen_glmm_data(
    mean.formula = ~ time*trt,
    re.formula = ~1 + time,
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
