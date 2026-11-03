source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
R <- readRDS('results_v4/S7d2/gates/GG3_rows.rds')
F <- R[is.finite(dense_stretches)]
F[, h := CI_WALK_FRAC * fifelse(is.finite(se_mod) & se_mod > 0, se_mod, 0.1)]
F[, c_ok := abs(dense_left_steps) * dense_step <= pmax(dense_step, h) & abs(dense_right_steps) * dense_step <= pmax(dense_step, h)]
B <- F[!(c_ok %in% TRUE)]
print(B[, .N, by = .(N, ntime)][order(N, ntime)])
TK <- s7d2_tasks()
res <- rbindlist(parallel::mclapply(seq_len(nrow(B)), mc.cores = 18, FUN = function(r) { b <- B[r]
  t <- unique(TK[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen); set.seed(FLIP_SEED(b$sim_id)); S <- flip_signs(P$N)
  w <- b$hi - b$lo; bs <- seq(b$lo, b$hi, length.out = 2001L)
  pv <- vapply(bs, function(x) { p <- flip_p_b(P, S, x)$p; if (is.na(p)) 1 else p }, 0)
  rj <- pv < 0.05
  data.table(N = b$N, T = b$ntime, est = b$est, lo = b$lo, hi = b$hi, h = b$h, share_rej_inside = mean(rj), n_islands = sum(diff(c(FALSE, rj)) == 1L),
             max_island_width = if (any(rj)) { r2 <- rle(rj); max(r2$lengths[r2$values]) * (bs[2] - bs[1]) } else 0, min_p_inside = min(pv),
             nref = nrow(S$G)) }))
print(res, digits = 3)
cat('share of the interval width that is rejected: median', median(res$share_rej_inside), ' max', max(res$share_rej_inside), '\n')
