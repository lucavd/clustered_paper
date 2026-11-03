# =============================================================================
# s7d2_gate_check.R -- gate GG3 (before the run) and the reported shape of the accepted set, on 1,080 datasets:
# 27 mechanism x scenario combinations (22 of the three study mechanisms, 5 of the copula), N 10 and 25, T 3 and 10
# (M1c: 3 and 6), replicates 1-5, H0 and H1.
# (a) at b0 = 0.1 and -0.3: reduced glm.fit with offset = glm(formula + offset(b0 * time * trt)) (coefficients within
#     1e-8), and the sum of effective score contributions = numDeriv gradient of the log-likelihood in the time:trt
#     coefficient with the other coefficients at the reduced fit (1e-6 relative or 1e-8 absolute);
# (b) at the centre |sum u| < 1e-8 and p = 1, in every dataset with an estimate;
# (c) (D-57) p on a grid of 2,001 values of b0 in [lo - width, hi + width] for finite intervals: passes if no rejected
#     run inside the reported interval is at least one walk step wide (the walk cannot miss such a run); reported: the
#     share of the interval that is rejected, the number of accepted stretches in the window, distance of the first
#     rejected value left / right of the centre from the found ends.
# Output: results_v4/S7d2/gates/GG3_rows.rds, tables/V4_S7d2_GG3.csv
# =============================================================================
source(path.expand("~/MMM/R/s7d2/s7d2_common.R")); suppressMessages(library(numDeriv))
TK <- s7d2_tasks(); u <- unique(TK[, .(dgp, scen, model)])
stopifnot(nrow(u) == 27L)
CB <- u[, .(N = c(10L, 25L)), by = .(dgp, scen, model)][, .(ntime = if (scen == "M1c") c(3L, 6L) else c(3L, 10L)), by = .(dgp, scen, model, N)]
CB <- CB[, .(label = c("mmlb", "manual")), by = .(dgp, scen, model, N, ntime)]
stopifnot(nrow(CB) == 216L)
FORM0 <- function(scen) switch(scen, M2 = binY ~ time + trt + age, M5 = binY ~ time + trt + dose, binY ~ time + trt)
one <- function(r) {
  t <- CB[r]; gen <- s7d2_gen(t); th <- truth_of(t$dgp, t$scen, t$ntime, t$label)
  rbindlist(lapply(1:5, function(i) {
    set_rep_seed(t$model, i); d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) NULL)
    if (is.null(d)) return(data.table(sim_id = i, generated = FALSE))
    d <- d[order(d$id, d$time), ]
    set.seed(FLIP_SEED(i)); ci <- flip_ci(d, t$scen, th)
    P <- flip_prep(d, t$scen); set.seed(FLIP_SEED(i)); S <- flip_signs(P$N)
    a <- rbindlist(lapply(c(0.1, -0.3), function(b) {
      r1 <- flip_p_b(P, S, b); if (is.null(r1$coef)) return(data.table(b = b, ok = FALSE))
      dd <- d; dd$off <- b * as.numeric(d$time) * as.numeric(d$trt)
      g <- tryCatch(suppressWarnings(glm(update(FORM0(t$scen), . ~ . + offset(off)), family = binomial, data = dd, control = GLMC)), error = function(e) NULL)
      dc <- if (is.null(g)) NA_real_ else max(abs(unname(coef(g)) - unname(r1$coef)))
      eta0 <- as.numeric(P$X0 %*% r1$coef)
      ll <- function(beta) { eta <- eta0 + beta * P$x; sum(P$y * eta - log1p(exp(eta))) }
      gr <- numDeriv::grad(ll, b)
      data.table(b = b, ok = TRUE, coef_diff = dc, su = r1$su, grad = gr)
    }))
    dense_n <- NA_integer_; dl <- NA_real_; dr <- NA_real_; run_in <- NA_real_; rej_in <- NA_real_
    if (is.finite(ci$lo) && is.finite(ci$hi)) {
      w <- ci$hi - ci$lo; bs <- seq(ci$lo - w, ci$hi + w, length.out = 2001L); stp <- bs[2] - bs[1]
      acc <- vapply(bs, function(b) !rejected(flip_p_b(P, S, b)), logical(1))
      dense_n <- sum(diff(c(FALSE, acc)) == 1L)
      L <- bs[bs < ci$est & !acc]; R <- bs[bs > ci$est & !acc]
      dl <- if (length(L)) (ci$lo - max(L)) / stp else NA_real_; dr <- if (length(R)) (min(R) - ci$hi) / stp else NA_real_
      ins <- bs > ci$lo & bs < ci$hi; rr <- rle(!acc[ins]); rej_in <- mean(!acc[ins])
      run_in <- if (any(rr$values)) (max(rr$lengths[rr$values]) - 1) * stp else 0   # lower bound of the widest rejected run inside
    }
    data.table(sim_id = i, generated = TRUE, se_mod = ci$se_mod, dense_step = if (is.finite(ci$lo) && is.finite(ci$hi)) 3 * (ci$hi - ci$lo) / 2000 else NA_real_, est = ci$est, est_why = ci$est_why, lo = ci$lo, hi = ci$hi, p_est = ci$p_est, su_est = ci$su_est,
               a_n = sum(a$ok), a_coef_max = suppressWarnings(max(a$coef_diff, na.rm = TRUE)),
               a_score_ok = all(with(a[ok == TRUE], abs(su - grad) <= pmax(1e-6 * abs(grad), 1e-8))),
               a_score_maxabs = suppressWarnings(max(abs(a$su - a$grad), na.rm = TRUE)),
               dense_stretches = dense_n, dense_left_steps = dl, dense_right_steps = dr,
               walk_h = walk_step(ci$se_mod), run_in = run_in, rej_in = rej_in)
  }))[, `:=`(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
}
W <- as.integer(Sys.getenv("GG3_WORKERS", "20"))
R <- rbindlist(parallel::mclapply(seq_len(nrow(CB)), one, mc.cores = W, mc.preschedule = FALSE), fill = TRUE)
dir.create(file.path(OUT_S7D2, "gates"), showWarnings = FALSE); saveRDS(R, file.path(OUT_S7D2, "gates/GG3_rows.rds"))
E <- R[generated == TRUE & is.finite(est) & is.na(est_why)]
F <- R[is.finite(dense_stretches)]
F[, c_ok := run_in < 1.01 * walk_h]
G <- data.table(
  part = c("datasets generated", "(a) offset fits compared", "(a) glm.fit vs glm, max coef diff", "(a) score = numDeriv, datasets failing",
           "(b) datasets with centre", "(b) |sum u| at centre >= 1e-8", "(b) p at centre != 1",
           "(c) finite intervals", "(c) rejected run inside the interval at least one walk step wide",
           "reported: intervals with a narrower rejected run inside", "reported: median share of the interval rejected, among those",
           "reported: max share of the interval rejected", "reported: >1 accepted stretch on the dense window", "reported: min p at centre"),
  value = c(R[generated == TRUE, .N], R[generated == TRUE, sum(a_n)], R[generated == TRUE, max(a_coef_max[is.finite(a_coef_max)])],
            R[generated == TRUE & a_n > 0, sum(a_score_ok %in% FALSE)], nrow(E), E[, sum(!(abs(su_est) < 1e-8), na.rm = TRUE) + sum(is.na(su_est))],
            E[, sum(is.na(p_est) | p_est != 1)], nrow(F), F[, sum(!(c_ok %in% TRUE))],
            F[, sum(rej_in > 0)], F[rej_in > 0, if (.N) median(rej_in) else 0], F[, max(rej_in)], F[, sum(dense_stretches > 1)], E[, min(p_est, na.rm = TRUE)]))
G[, pass := c(NA, NA, value[3] < 1e-8, value[4] == 0, NA, value[6] == 0, value[7] == 0, NA, value[9] == 0, NA, NA, NA, NA, NA)]
print(G); fwrite(G, file.path(TABDIR, "V4_S7d2_GG3.csv"))
cat("GG3", if (all(G$pass, na.rm = TRUE)) "PASS" else "FAIL", "\n")
