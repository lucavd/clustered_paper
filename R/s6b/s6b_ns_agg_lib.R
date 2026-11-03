# =============================================================================
# s6b_ns_agg_lib.R -- S6b part B: aggregation with the definitions of R/s5/s5_build.R, R/s6/s6_dp_rates.R and
# R/s6/s6_rule_safety.R (own set, usability U0 and not a degenerate design; normal and t reference with N - 4
# degrees of freedom; Bradley limits with the binomial MCSE of the design point). Checked against the study tables
# on M3 (gate GB4) before it is applied to the new scenarios.
# =============================================================================
NP3 <- 4L
conv_rows <- function(R) {     # chunk rows -> rows_v4 layout (as R/s4/s4_load.R)
  R <- R[method != "gee_ss"]
  B <- R[method == LRT & !is.na(term)]
  B[, `:=`(method = SENS, p = 2 * pnorm(-abs(estimate / se)), stat = estimate / se)]
  R <- rbindlist(list(R, B), use.names = TRUE, fill = TRUE); add_usability(R); R
}
# V: dgp scen N ntime label method sim_id estimate se p U0 deg th tok ; ND: dgp scen N ntime nd
indicators <- function(V) {
  V[, u0 := U0 & !deg]; V[, arm := fifelse(label == "mmlb", "H1", "H0")]
  V[, r := u0 & fifelse(is.na(p), FALSE, p < 0.05)]
  V[, rt := u0 & method != LRT & fifelse(is.finite(estimate / se), 2 * pt(-abs(estimate / se), df = N - NP3) < 0.05, FALSE)]
  V[, cv := u0 & tok & fifelse(is.finite(th), abs(estimate - th) <= Z975 * se, FALSE)]
  V[, er := fifelse(u0 & tok, estimate - th, NA_real_)]
  invisible(V)
}
dp_rates <- function(V, ND) {
  X <- V[, .(n = sum(u0), k = sum(r), kt = sum(rt), kc = sum(cv), bias = mean(er[u0], na.rm = TRUE), empse = sd(estimate[u0]),
             modse = mean(se[u0]), modse_rms = sqrt(mean(se[u0]^2))), by = .(dgp, scen, N, ntime, arm, method)]
  # full layout, as the META table of R/s5/s5_build.R: a method with no row at a design point has n = 0
  G <- ND[, .(dgp, scen, N, ntime, nd)][rep(seq_len(.N), each = 2L * uniqueN(V$method))]
  G[, `:=`(arm = rep(rep(c("H0", "H1"), each = uniqueN(V$method)), .N / (2L * uniqueN(V$method))), method = rep(sort(unique(V$method)), .N / uniqueN(V$method)))]
  X <- merge(G, X, by = c("dgp", "scen", "N", "ntime", "arm", "method"), all.x = TRUE)
  for (v in c("n", "k", "kt", "kc")) set(X, which(is.na(X[[v]])), v, 0)
  X[method == LRT, kt := NA_real_]; X
}
rule_safety_dp <- function(X) {
  L <- rbind(X[, .(dgp, scen, N, ntime, arm, method, n, nd, kk = k, ref = "n")],
             X[method != LRT, .(dgp, scen, N, ntime, arm, method, n, nd, kk = kt, ref = "t")])
  L[, strategy := paste(method, ref, sep = "|")]
  H0 <- L[arm == "H0"][, rate := fifelse(n > 0, 100 * kk / n, NA_real_)]
  H0[, mc := fifelse(n > 0, 100 * sqrt((rate / 100) * (1 - rate / 100) / n), NA_real_)]
  H0[, `:=`(assess = n >= 20, pu = 100 * n / nd)]
  H0[, `:=`(hi = assess & rate - 2 * mc > 7.5, lo = assess & rate + 2 * mc < 2.5)]
  H1 <- L[arm == "H1", .(dgp, scen, N, ntime, strategy, upow = 100 * kk / nd)]
  D <- merge(H0, H1, by = c("dgp", "scen", "N", "ntime", "strategy"))
  D[, band := fifelse(ntime <= 3, "T2-3", fifelse(ntime <= 6, "T4-6", "T8-10"))]
  D[, .(dgp, scen, N, ntime, band, strategy, n, nd, rate, mc, assess, hi, lo, pu, upow)]
}
# block bootstrap of means over cells (as R/s6/s6cc_compare.R): D has dgp scen model N ntime label sim_id + 0/1 columns
boot_counts <- function(M) { Cm <- matrix(0, 1000L, nrow(M)); for (b in seq_len(nrow(M))) Cm[, b] <- tabulate(M[b, ], 1000L); Cm }
cellrates <- function(D, ucol, rcol, CNT) {
  out <- list()
  for (m in unique(D$model)) {
    d <- D[model == m]; d[, cell := paste(dgp, scen, N, ntime, label)]; cl <- unique(d$cell)
    U <- matrix(0, length(cl), 1000L); R <- U; ix <- cbind(match(d$cell, cl), d$sim_id)
    U[ix] <- as.numeric(d[[ucol]]); R[ix] <- as.numeric(d[[rcol]])
    Wm <- cbind(1, CNT[[as.character(m)]]); out[[as.character(m)]] <- list(cell = cl, rate = 100 * (R %*% Wm) / (U %*% Wm))
  }
  list(cell = unlist(lapply(out, `[[`, "cell")), rate = do.call(rbind, lapply(out, `[[`, "rate")))
}
meanboot <- function(rt) { m <- colMeans(rt, na.rm = TRUE); c(est = m[1], mcse = sd(m[-1])) }
judge <- function(est, mc, lo, hi) { a <- est - 2 * mc; b <- est + 2 * mc
  fifelse(a >= lo & b <= hi, "confermata", fifelse(b < lo | a > hi, "smentita", "indeterminata")) }
# GLMM-native marginal truth: exact Gauss-Hermite quadrature on the weighted design (method of R/s2/s2_c_truth.R, M3 design)
gh_std <- function(q) { g <- fastGHQuad::gaussHermiteData(q); list(x = sqrt(2) * g$x, w = g$w / sqrt(pi)) }
marg_prob <- function(eta, v, gh) vapply(seq_along(eta), function(k) sum(gh$w * plogis(eta[k] + sqrt(v[k]) * gh$x)), 0)
irls <- function(X, p, w) { beta <- rep(0, ncol(X))
  for (it in 1:500) { mu <- plogis(drop(X %*% beta)); W <- w * mu * (1 - mu)
    step <- solve(crossprod(X, X * W), crossprod(X, w * (p - mu))); beta <- beta + drop(step); if (max(abs(step)) < 1e-14) break }
  setNames(beta, colnames(X)) }
native_truth3 <- function(bc, D, T, q = 80) {
  gh <- gh_std(q); g <- CJ(time = 0:(T - 1), trt = 0:1)
  g[, w := fifelse(time == 0, 1, 0.7 + 0.3 * (T - 1 - time) / (T - 1)) * 0.5]
  X <- cbind(`(Intercept)` = 1, time = g$time, trt = g$trt, `time:trt` = g$time * g$trt); Z <- cbind(1, g$time)
  v <- rowSums((Z %*% D) * Z); unname(irls(X, marg_prob(drop(X %*% bc), v, gh), g$w)["time:trt"])
}
