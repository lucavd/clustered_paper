# V4 (S3, D-13/D-14): scenario M1c for GLMM-native. Design code identical to
# manual/generation_model1.R; beta_c = Zeger factor x marginal vector (M1 convention).
generate_data_model1c <- function(N, ntime, sim_id, methods = c("mmlb", "manual")) {
  methods <- match.arg(methods)
  beta_c <- M1C_ZEGER_F * m1c_beta_marginal(ntime, methods)

  nclust <- ifelse(
    runif(N) < 0.7,
    ntime,
    sample(seq_len(ntime - 1), N, replace = TRUE)
  )
  id <- rep(seq(N), nclust)
  trt <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  data <- data.frame(id, time, trt) |> arrange(id, time)

  gen_glmm_data(
    mean.formula = ~ time*trt,
    re.formula = ~1,
    beta_c = beta_c,
    D = M1C_VAR,
    id = id,
    data = data,
    Yname = "binY"
  ) |>
    dplyr::mutate(sim_id = sim_id)
}
