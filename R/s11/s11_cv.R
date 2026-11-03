# =============================================================================
# s11_cv.R -- D-90 helper (S11, D-91): in the convergence check of the marginalized model, a fit at the boundary is
# replaced by the fit without the slope terms only if nlm reports convergence for it (code 1 or 2) and its estimate and
# SE are finite; otherwise the dataset gives no usable result for the MMM with the check. The keys (dgp, scen, N, ntime,
# label, sim_id) of these datasets are built once by R/s11/s11_cv_keys.R from the stored check fits and saved in
# results_v4/S11/cv_keys.rds. s11_drop_cv(X, col, rows) sets X[[col]] to FALSE on those keys, ONLY on the rows selected
# by `rows` (the strategy with the check), by reference, and returns how many TRUE values it switched off.
# Every call is marked "# D-90" in the calling script.
# =============================================================================
S11_KEYS_FILE <- "/home/user/MMM/results_v4/S11/cv_keys.rds"
S11_K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
s11_cv_keys <- local({ cache <- NULL; function() {
  if (is.null(cache)) { k <- data.table::as.data.table(readRDS(S11_KEYS_FILE)); stopifnot(identical(names(k)[1:6], S11_K6)); cache <<- k[, S11_K6, with = FALSE] }
  cache } })
s11_drop_cv <- function(X, col = "use", rows, keys = s11_cv_keys(), what = deparse(substitute(X))) {
  stopifnot(data.table::is.data.table(X), all(S11_K6 %in% names(X)), col %in% names(X), is.logical(rows), length(rows) == nrow(X))
  kk <- data.table::copy(keys)
  for (v in c("N", "ntime", "sim_id")) data.table::set(kk, j = v, value = methods::as(kk[[v]], class(X[[v]])[1]))
  ix <- X[kk, on = S11_K6, which = TRUE, nomatch = NULL]
  ix <- ix[rows[ix] %in% TRUE]
  n_off <- sum(X[[col]][ix] %in% TRUE)
  if (length(ix)) data.table::set(X, i = ix, j = col, value = FALSE)
  message(sprintf("D-90: %s$%s switched off on %d of %d rows of the MMM with the check on a no-result key", what, col, n_off, length(ix)))
  invisible(n_off)
}
