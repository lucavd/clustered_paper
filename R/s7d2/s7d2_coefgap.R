# s7d2_coefgap.R -- distribution of max |coefficient| of the reduced fit with offset over the dense grids of the GG3
# datasets with a finite interval (cold starts), to place the divergence threshold of D-57 in the gap.
source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
R <- readRDS('results_v4/S7d2/gates/GG3_rows.rds'); F <- R[is.finite(dense_stretches)]
TK <- s7d2_tasks()
res <- rbindlist(parallel::mclapply(seq_len(nrow(F)), mc.cores = 20, FUN = function(r) { b <- F[r]
  t <- unique(TK[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen)
  bs <- seq(b$lo - (b$hi - b$lo), b$hi + (b$hi - b$lo), length.out = 401L)
  mc <- vapply(bs, function(x) { f <- tryCatch(suppressWarnings(glm.fit(P$X0, P$y, family = binomial(), offset = x * P$x, control = GLMC)), error = function(e) NULL)
                                 if (is.null(f) || !all(is.finite(f$coefficients))) NA_real_ else max(abs(f$coefficients)) }, 0)
  data.table(N = b$N, T = b$ntime, mc = mc) }))
res[, bin := cut(log10(mc), c(-Inf, 0, 1, 1.5, 2, 3, 4, 6, 10, 14, Inf))]
print(res[, .N, by = bin][order(bin)])
print(res[, .(share_ge_1e3 = mean(mc >= 1e3, na.rm = TRUE), max_below_1e3 = max(mc[mc < 1e3], na.rm = TRUE), min_above_1e3 = min(mc[mc >= 1e3], na.rm = TRUE)), by = .(N)])
