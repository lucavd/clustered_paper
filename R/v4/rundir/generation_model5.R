# V4 (S3, D-01): copy of .rundir/generation_model5.R with the H0 sigma/gamma set to the H1 values.
generate_data_model5 <- function(
    N,
    ntime,
    sim_id,
    methods = c("mmlb", "manual")
) {
  methods <- match.arg(methods)
  # H0 branch differs from H1 in the tested parameter only (see generation_model1.R).
  beta_vec <- switch(
    methods,
    # (Intercept), time, trt, dose, time:trt -- interaction LAST (R term order)
    mmlb   = c(-2, 0.1, 0.2, 0.05, 0.2),
    manual = c(-2, 0.1, 0.2, 0.05, 0)
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
  sex = rep(rbinom(N, size = 1, prob = 0.5), nclust) # Binary sex variable
  steroids = rep(sample(1:3, N, replace = TRUE), nclust) # Steroid treatment levels
  # dose VARIES within subject -- unlike age and sex, which are drawn once per
  # subject. A random slope needs a covariate that varies within the cluster,
  # otherwise its column is proportional to the intercept's and the two are not
  # separately identified.
  dose <- as.numeric(scale(rnorm(length(id))))
  data <- data.frame(id, time, trt, age, sex, steroids, dose) |> arrange(id, time)

   beta_vec <- beta_vec

   # lv.formula = ~time + dose has three columns, so sigma has three elements
   sigma <- switch(
    methods,
    mmlb = c(1, 0.01, 0.01),
    manual = c(1, 0.01, 0.01))

   gamma <- switch(
     methods,
     mmlb = c(1, 0.1, 0.05),
     manual = c(1, 0.1, 0.05))

   GenBinaryY(
    mean.formula = ~ time*trt + dose,
    lv.formula = ~time + dose,
    t.formula = ~trt + steroids,
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
