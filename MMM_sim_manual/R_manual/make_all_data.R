make_all_data <- function(n, ntime, sim_id, model, method) {
  generate_data(
    model  = model,
    method = method,
    N      = n,
    ntime  = ntime,
    sim_id   = sim_id
  )
}
