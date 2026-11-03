generate_data_model1 <- function(
    N,
    ntime,
    sim_id,
    methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # H0 branch differs from H1 in the TESTED parameter only; intercept and main
  # effects are held fixed so type I error and power are evaluated at the same
  # outcome prevalence (previously c(0,0,0,0) also moved the intercept). 2026-08-05.
  beta_c <- switch(
    methods,
    # (Intercept), time, trt, time:trt -- conditional; marginal ≈ (-2, 0.1, 0.2, 0.2)
    mmlb   = c(-2.32, 0.116, 0.232, 0.232),
    manual = c(-2.32, 0.116, 0.232, 0)
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

  D <- 1.0  # random intercept variance

  gen_glmm_data(
    mean.formula = ~ time*trt,
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
