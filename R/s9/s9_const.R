# =============================================================================
# s9_const.R -- D-82 helper (S9, D-83): a non-degenerate dataset whose outcome is constant (sum_y == 0 or
# sum_y == n_obs) gives no usable result for any method. The keys (dgp, scen, N, ntime, label, sim_id) of these
# datasets are built once by R/s9/s9_const_keys.R from the recs of the S7d run (three study mechanisms, same seeds
# as the stored study rows) and of the S7e run (Copula-AR1), and saved in results_v4/S9/const_keys.rds.
# s9_drop_const(X, col) sets X[[col]] to FALSE on those keys (by reference) and returns how many TRUE values it
# switched off. Every call is marked "# D-82" in the calling script.
# =============================================================================
S9_KEYS_FILE <- "/home/user/MMM/results_v4/S9/const_keys.rds"
S9_K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
s9_const_keys <- local({ cache <- NULL; function() {
  if (is.null(cache)) { k <- data.table::as.data.table(readRDS(S9_KEYS_FILE)); stopifnot(identical(names(k), S9_K6)); cache <<- k }
  cache } })
s9_drop_const <- function(X, col = "use", keys = s9_const_keys(), what = deparse(substitute(X))) {
  stopifnot(data.table::is.data.table(X), all(S9_K6 %in% names(X)), col %in% names(X))
  kk <- data.table::copy(keys)
  for (v in c("N", "ntime", "sim_id")) data.table::set(kk, j = v, value = methods::as(kk[[v]], class(X[[v]])[1]))
  ix <- X[kk, on = S9_K6, which = TRUE, nomatch = NULL]
  n_off <- sum(X[[col]][ix] %in% TRUE)
  if (length(ix)) data.table::set(X, i = ix, j = col, value = FALSE)
  message(sprintf("D-82: %s$%s switched off on %d of %d constant-outcome rows", what, col, n_off, length(ix)))
  invisible(n_off)
}
