# s7d2_diag17.R -- diagnosis of the GG3 (c) failures with the fixed-step walk: same datasets, walk recomputed with
# warm starts (as in s7d2_ci.R) and with cold starts at every b0; p at the widest rejected run inside, cold vs warm.
source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
R <- readRDS('results_v4/S7d2/gates/GG3_rows.rds'); F <- R[is.finite(dense_stretches)]
B <- F[!(run_in < 1.01 * walk_h)]
TK <- s7d2_tasks()
res <- rbindlist(parallel::mclapply(seq_len(nrow(B)), mc.cores = 17, FUN = function(r) { b <- B[r]
  t <- unique(TK[dgp == b$dgp & scen == b$scen & N == b$N & ntime == b$ntime & label == b$label, .(dgp, scen, N, ntime, label, model)])
  gen <- s7d2_gen(t); set_rep_seed(t$model, b$sim_id); d <- as.data.frame(gen(t$N, t$ntime, b$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  P <- flip_prep(d, t$scen); set.seed(FLIP_SEED(b$sim_id)); S <- flip_signs(P$N)
  bs <- seq(b$lo, b$hi, length.out = 2001L)
  cold <- lapply(bs, function(x) flip_p_b(P, S, x))
  rj <- vapply(cold, rejected, TRUE); rr <- rle(rj); k <- which.max(ifelse(rr$values, rr$lengths, 0)); end <- cumsum(rr$lengths)[k]; st <- end - rr$lengths[k] + 1
  mid <- bs[floor((st + end) / 2)]
  # walk points that fall inside the widest run
  h <- b$walk_h; kk <- seq_len(ceiling(20 / h)); wp <- c(b$est - kk * h, b$est + kk * h); wp <- wp[wp >= bs[st] & wp <= bs[end]]
  pw <- if (length(wp)) vapply(wp, function(x) flip_p_b(P, S, x)$p, 0) else NA_real_
  # p at the run midpoint, cold, and warm from the coefficients of the centre
  cen <- flip_p_b(P, S, b$est)
  warm <- flip_p_b(P, S, mid, cen$coef); cm <- flip_p_b(P, S, mid)
  data.table(dgp = b$dgp, scen = b$scen, N = b$N, T = b$ntime, label = b$label, sim = b$sim_id, est = b$est, lo = b$lo, hi = b$hi, h = h,
             run_from = bs[st], run_to = bs[end], n_walk_in_run = length(wp), p_walk_pts_cold = paste(signif(head(pw, 4), 3), collapse = ' '),
             p_mid_cold = cm$p, p_mid_warm = warm$p, coef_diff = if (!is.null(warm$coef) && !is.null(cm$coef)) max(abs(warm$coef - cm$coef)) else NA_real_,
             maxcoef_cold = if (!is.null(cm$coef)) max(abs(cm$coef)) else NA_real_) }))
print(res, digits = 3)
