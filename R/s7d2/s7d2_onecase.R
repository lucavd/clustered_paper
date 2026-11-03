# s7d2_onecase.R -- one dry-run dataset (GLMM-native M6, N 15, T 4, H1, replicate 70) whose control arm has only 0 outcomes: reduced fit, fitted range, sum of contributions and p at several b0 (D-57).
source(path.expand('~/MMM/R/s7d2/s7d2_common.R'))
TK <- s7d2_tasks(); t <- unique(TK[dgp == 'GLMM-native' & scen == 'M6' & N == 15 & ntime == 4 & label == 'mmlb', .(dgp, scen, N, ntime, label, model)])
gen <- s7d2_gen(t); i <- 70L; set_rep_seed(t$model, i); d <- as.data.frame(gen(15, 4, i, 'mmlb')); d <- d[order(d$id, d$time), ]
P <- flip_prep(d, 'M6'); set.seed(FLIP_SEED(i)); S <- flip_signs(P$N)
cat('n obs', length(P$y), 'sum y', sum(P$y), ' y among x>0:', sum(P$y[P$x > 0]), '/', sum(P$x > 0), '\n')
print(table(trt = d$trt, time = d$time, y = d$binY))
for (b in c(-1.5, -1.8, -1.8253, -1.82531, -1.8254, -1.83, -2, -3)) {
  f <- suppressWarnings(glm.fit(P$X0, P$y, family = binomial(), offset = b * P$x, control = GLMC))
  r <- flip_p_b(P, S, b)
  mu <- f$fitted.values; e <- P$y - mu
  cat(sprintf('b %.5f conv %s iter %d coef %s | mu range %.2e %.2e | sum u %.3e p %.4f\n', b, f$converged, f$iter, paste(signif(f$coefficients, 4), collapse = ','), min(mu), max(mu), r$su, r$p))
}
