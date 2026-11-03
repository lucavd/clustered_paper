# s7d2_zerowidth.R -- finite intervals of (near) zero width: counts by N and T, and the p-values around the centre
# for the first cases (centre, centre +/- one walk step).
source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
TA <- s7d2_tasks()
Z <- rbindlist(parallel::mclapply(TA$file, function(f) { x <- readRDS(f); x[is.finite(lo) & is.finite(hi) & (hi - lo) < 1e-3] }, mc.cores = 16L), fill = TRUE)
cat('finite intervals with width < 1e-3:', nrow(Z), ' width exactly 0:', Z[hi == lo, .N], '\n')
print(Z[, .N, by = .(N, ntime)][order(N, ntime)])
print(summary(Z$p_est))
print(Z[, .N, by = .(p_est_na = is.na(p_est), p_est_is1 = p_est %in% 1)])
for (r in 1:3) { b <- Z[r]; t <- unique(TA[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen); set.seed(FLIP_SEED(b$sim_id)); S <- flip_signs(P$N); h <- walk_step(b$se_mod)
  pp <- sapply(b$est + c(-2, -1, -0.5, -0.01, 0, 0.01, 0.5, 1, 2) * h, function(x) { r <- flip_p_b(P, S, x); if (is.na(r$p)) r$why else sprintf('%.4f', r$p) })
  cat(sprintf('%s %s N%d T%d %s i%d est %.4f se %.3g lo %.6f hi %.6f | p at est + (-2,-1,-.5,-.01,0,.01,.5,1,2) h: %s\n', b$dgp, b$scen, b$N, b$ntime, b$label, b$sim_id, b$est, b$se_mod, b$lo, b$hi, paste(pp, collapse = ' ')))
  cat('  y by arm:', tapply(P$y, d$trt, sum), ' n by arm:', tapply(P$y, d$trt, length), '\n') }
