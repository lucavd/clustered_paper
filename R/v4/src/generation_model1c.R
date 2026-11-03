# V4 (S3, D-13/D-14): scenario M1c for GenBinaryY. Design code identical to
# generation_model1.R (same random-number sequence); parameters from m1c_params.R.
generate_data_model1c <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  beta_vec <- m1c_beta_marginal(ntime, methods)

  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id <- rep(seq(N), nclust)
  trt <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  data <- data.frame(id, time, trt) |> arrange(id, time)

  GenBinaryY(
    mean.formula = ~ time*trt,
    lv.formula = ~1,
    t.formula = ~1,
    beta = beta_vec,
    sigma = M1C_SIGMA_GBY,
    gamma = M1C_GAMMA_GBY,
    id = id,
    data = data,
    q = 20,
    Yname = "binY"
  ) |>
    dplyr::mutate(sim_id = sim_id)
}
