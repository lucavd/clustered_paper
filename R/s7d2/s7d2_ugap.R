# s7d2_ugap.R -- distribution of max_i |u_i| (effective score contributions) of the reduced fit with offset, over 401
# values of b0 per GG3 dataset with an estimate (window: estimate +/- 5), to place the "score numerically zero" rule of
# D-57 in the gap.
source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
R <- readRDS('results_v4/S7d2/gates/GG3_rows.rds'); F <- R[generated == TRUE & is.finite(est) & is.na(est_why)]
TK <- s7d2_tasks()
res <- rbindlist(parallel::mclapply(seq_len(nrow(F)), mc.cores = 20, FUN = function(r) { b <- F[r]
  t <- unique(TK[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen)
  mu_ <- vapply(seq(b$est - 5, b$est + 5, length.out = 401L), function(x) {
    f <- tryCatch(suppressWarnings(glm.fit(P$X0, P$y, family = binomial(), offset = x * P$x, control = GLMC)), error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$converged) || !all(is.finite(f$coefficients))) return(NA_real_)
    mu <- f$fitted.values; e <- P$y - mu; w <- mu * (1 - mu)
    proj <- tryCatch(solve(crossprod(P$X0, P$X0 * w), colSums(P$X0 * (w * P$x))), error = function(e) NULL); if (is.null(proj)) return(NA_real_)
    max(abs(rowsum(P$x * e, P$id, reorder = TRUE)[, 1] - rowsum(P$X0 * e, P$id, reorder = TRUE) %*% proj)) }, 0)
  data.table(k = r, N = b$N, T = b$ntime, mu = mu_) }))
res[, bin := cut(log10(mu), c(-Inf, -14, -12, -10, -8, -6, -4, -2, 0, Inf))]
print(res[, .N, by = bin][order(bin)])
cat('datasets with any value < 1e-8:', res[, any(mu < 1e-8, na.rm = TRUE), by = k][, sum(V1)], ' of', uniqueN(res$k), '\n')
print(res[mu < 1e-8, .N, by = .(N, T)])
