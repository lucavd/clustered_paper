# s7d2_fuzzy.R -- dry-run cases where the point 1e-4 beyond a reported end is accepted: p on 401 values within +/- 1e-3
# of the end (cold starts), number of crossings of 0.05 and the p-values met.
source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
X <- rbindlist(lapply(list.files('results_v4/S7d2/rows', full.names = TRUE), readRDS), fill = TRUE)
B <- X[is.finite(est) & is.na(est_why) & ((is.finite(lo) & !(p_lo_out < 0.05)) | (is.finite(hi) & !(p_hi_out < 0.05)))]
TK <- s7d2_tasks()
for (r in seq_len(nrow(B))) { b <- B[r]
  t <- unique(TK[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen); set.seed(FLIP_SEED(b$sim_id)); S <- flip_signs(P$N)
  e <- if (is.finite(b$lo) && !(b$p_lo_out < 0.05)) b$lo else b$hi; dir <- if (identical(e, b$lo)) -1 else 1
  bs <- e + dir * seq(-1e-3, 1e-3, length.out = 401L); pv <- vapply(bs, function(x) flip_p_b(P, S, x)$p, 0)
  cr <- sum(diff(pv < 0.05) != 0)
  cat(sprintf('%s %s N%d T%d %s i%d end %.5f side %+d: crossings of 0.05 within +/-1e-3 = %d; p values met: %s\n', b$dgp, b$scen, b$N, b$ntime, b$label, b$sim_id, e, dir, cr,
              paste(sort(unique(signif(pv, 3))), collapse = ' ')))
}
