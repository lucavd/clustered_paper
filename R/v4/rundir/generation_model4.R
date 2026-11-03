# V4 (S3, D-01): copy of .rundir/generation_model4.R with the H0 sigma/gamma set to the H1 values.
generate_data_model4 <- function(
    N,
    ntime,
    sim_id,
    methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # M4 deliberately uses a small effect (time:trt = 0.05) to probe a
  # near-null regime. H0 branch differs in the tested parameter only.
  beta_vec <- switch(
    methods,
    # (Intercept), time, trt, time:trt
    mmlb   = c(-2, 0.2, 0.2, 0.05),
    manual = c(-2, 0.2, 0.2, 0)
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
    mmlb = c(1, 0.01),
    manual = c(1, 0.01) ) # D-01: H0 sigma = H1 sigma

  gamma <- switch(
    methods,
    mmlb = c(1, 0.01),
    manual = c(1, 0.01) )

  GenBinaryY(
    mean.formula = ~ time*trt,
    lv.formula = ~time,
    t.formula = ~trt,
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
