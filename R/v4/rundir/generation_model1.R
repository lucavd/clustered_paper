generate_data_model1 <- function(
    N,
    ntime,
    sim_id,
    methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # H0 branch ("manual") differs from H1 ("mmlb") in the tested parameter ONLY.
  # It previously used c(0,0,0,0), which also moved the intercept from -2 to 0 and
  # thus evaluated type I error at ~50% prevalence while power was evaluated at
  # ~12%. Intercept and main effects are now held fixed. Fixed 2026-08-05.
  beta_vec <- switch(
    methods,
    # (Intercept), time, trt, time:trt
    mmlb   = c(-2, 0.1, 0.2, 0.2),
    manual = c(-2, 0.1, 0.2, 0)
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
  data <- data.frame(id, time, trt) |> arrange(id, time)

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
    mean.formula = ~ time*trt,
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
