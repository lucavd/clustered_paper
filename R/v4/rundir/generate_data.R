generate_data <- function(model, method, N, ntime, sim_id) {
  model |>
    switch(
      "1" = generate_data_model1(N, ntime, sim_id, method),
      "2" = generate_data_model2(N, ntime, sim_id, method),
      "3" = generate_data_model3(N, ntime, sim_id, method),
      "4" = generate_data_model4(N, ntime, sim_id, method),
      "5" = generate_data_model5(N, ntime, sim_id, method),
      stop("Model or method not recognized")
    )
}
