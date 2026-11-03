# =============================================================================
# s7d2_ci.R -- S7d2: confidence interval for time:trt by inverting the sign-flip score test (D-52; design fixed in
# docs/steps/S7d2.md before any calculation). Test at b0 = the primary test of R/s7d/s7d_flip.R with the null
# time:trt = b0 imposed by the offset b0 * time * trt in the reduced model (glm.fit, epsilon 1e-12, maxit 100).
# One sign matrix per dataset, the same for every b0, drawn exactly as test_flip() draws it after
# set.seed(20261022 + i): all 2^N signs if N <= 12, else 1,999 random rows (+ identity in the p-value).
# Centre = independence GEE estimate (full logistic glm.fit), where the effective score sums to 0 and p = 1.
# Search per side: see ci_side() below (walk of D-57); not rejected at the bound -> that side is unbounded. A b0 at which the
# reduced fit fails, diverges (max |coefficient| >= 1e3) or has a numerically zero score (max |u_i| < 1e-6) (D-57)
# counts as "not rejected" (conservative) and is counted. Neither rule applies to the reproduction at b0 = 0.
# =============================================================================
CI_ALPHA <- 0.05; CI_BMAX <- 10; CI_TOL <- 1e-4; CI_EPS_BEYOND <- 1e-4
CI_DIVERGE <- 1e3   # D-57: reduced or full fit with max |coefficient| >= 1e3 has diverged (separation): not evaluable
CI_UZERO <- 1e-6    # D-57: max_i |u_i| < 1e-6 -> score numerically zero (time:trt not identified, e.g. an arm all 0): not evaluable
GLMC <- glm.control(epsilon = 1e-12, maxit = 100)

flip_prep <- function(d, scen) {   # d sorted by id, time
  list(id = as.integer(factor(d$id)), y = as.numeric(d$binY), x = as.numeric(d$time) * as.numeric(d$trt),
       X0 = s7d_reduced_X(d, scen), N = length(unique(d$id)))
}
flip_signs <- function(N) {        # RNG state set by the caller; identical draw to test_flip()
  if (N <= ENUM_FLIP_MAXN) list(G = all_signs(N), enum = TRUE)
  else list(G = matrix(sample(c(-1, 1), B_FLIP * N, replace = TRUE), B_FLIP, N), enum = FALSE)
}
# p-value of the sign-flip test of time:trt = b at fixed signs S; returns p, reduced coefficients, sum of u
flip_p_b <- function(P, S, b, start = NULL, diverge = TRUE) {   # diverge = FALSE only for the reproduction at b0 = 0
  f <- tryCatch(suppressWarnings(glm.fit(P$X0, P$y, family = binomial(), offset = b * P$x, start = start, control = GLMC)),
                error = function(e) NULL)
  ok <- !is.null(f) && isTRUE(f$converged) && all(is.finite(f$coefficients)) && all(is.finite(f$fitted.values))
  dv <- ok && diverge && max(abs(f$coefficients)) >= CI_DIVERGE
  if (!ok || dv) { if (!is.null(start)) return(flip_p_b(P, S, b, NULL, diverge))
    return(list(p = NA_real_, coef = NULL, su = NA_real_, why = if (dv) "reduced fit diverged" else "reduced fit failed")) }
  mu <- f$fitted.values; e <- P$y - mu; w <- mu * (1 - mu)
  I00 <- crossprod(P$X0, P$X0 * w); Ix0 <- colSums(P$X0 * (w * P$x))
  proj <- tryCatch(solve(I00, Ix0), error = function(e) NULL)
  if (is.null(proj)) return(list(p = NA_real_, coef = f$coefficients, su = NA_real_, why = "information singular"))
  u <- as.numeric(rowsum(P$x * e, P$id, reorder = TRUE)[, 1] - rowsum(P$X0 * e, P$id, reorder = TRUE) %*% proj)
  if (!(sum(u^2) > 0)) return(list(p = NA_real_, coef = f$coefficients, su = sum(u), why = "zero score variance"))
  if (diverge && max(abs(u)) < CI_UZERO) return(list(p = NA_real_, coef = f$coefficients, su = sum(u), why = "score numerically zero"))
  To <- abs(sum(u)); Ts <- abs(as.numeric(S$G %*% u))   # same arithmetic as test_flip()
  p <- if (S$enum) mean(Ts >= To * (1 - 1e-10)) else (1 + sum(Ts >= To * (1 - 1e-10))) / (B_FLIP + 1)
  list(p = p, coef = f$coefficients, su = sum(u), why = NA_character_)
}
rejected <- function(r) !is.na(r$p) && r$p < CI_ALPHA
walk_step <- function(s) CI_WALK_FRAC * min(if (is.finite(s) && s > 0) s else 0.1, CI_S_MAX)
# Search (D-57, correction before the production run): walk outward from the centre in fixed steps of
# h = min(s, 1) / 10 (s = model-based SE of the full fit, 0.1 if not finite) up to |b0| = 10; at the first rejected
# b0, bisection between it and the last accepted value until the gap is < 1e-4. The reported end is the accepted side:
# the stretch of accepted values around the estimate, up to rejected runs narrower than h. (The pre-registered
# doubling s, 2s, 4s, ... jumped over rejected regions into outer accepted regions: reports/morris/S7d2_logs/
# GG3_doubling.txt; a walk whose step grew 5% after 100 steps did the same far from the centre: GG3_walkgrow.txt.)
CI_WALK_FRAC <- 0.1; CI_S_MAX <- 1
ci_side <- function(P, S, bh, s, dir) {   # every fit from glm.fit's default start (D-57): p is a function of b0 only
  a <- bh; ne <- 0L; nna <- 0L; rej <- NA_real_; h <- walk_step(s)
  for (k in 1:(ceiling(2 * CI_BMAX / h) + 2)) {
    b <- bh + dir * k * h; atb <- FALSE
    if (dir > 0 && b >= CI_BMAX) { b <- CI_BMAX; atb <- TRUE }
    if (dir < 0 && b <= -CI_BMAX) { b <- -CI_BMAX; atb <- TRUE }
    r <- flip_p_b(P, S, b); ne <- ne + 1L
    if (is.na(r$p)) nna <- nna + 1L
    if (rejected(r)) { rej <- b; break }
    a <- b
    if (atb) return(list(end = dir * Inf, rej = NA_real_, ne = ne, nna = nna))
  }
  if (is.na(rej)) return(list(end = dir * Inf, rej = NA_real_, ne = ne, nna = nna))
  while (abs(rej - a) > CI_TOL) {
    m <- (a + rej) / 2; r <- flip_p_b(P, S, m); ne <- ne + 1L
    if (is.na(r$p)) nna <- nna + 1L
    if (rejected(r)) rej <- m else a <- m
  }
  list(end = a, rej = rej, ne = ne, nna = nna)
}
# full interval record for one dataset; RNG state must be set (set.seed(FLIP_SEED(i))) by the caller
flip_ci <- function(d, scen, theta) {
  P <- flip_prep(d, scen); S <- flip_signs(P$N)
  out <- list(p0_new = NA_real_, p0_guard = NA_real_, why0_guard = NA_character_, p_theta = NA_real_, theta_why = NA_character_, est = NA_real_, se_mod = NA_real_, est_why = NA_character_,
              lo = NA_real_, hi = NA_real_, rej_lo = NA_real_, rej_hi = NA_real_, p_lo = NA_real_, p_hi = NA_real_,
              p_lo_out = NA_real_, p_hi_out = NA_real_, p_est = NA_real_, su_est = NA_real_, n_eval = 0L, n_na = 0L)
  out$p0_new <- flip_p_b(P, S, 0, diverge = FALSE)$p
  r0 <- flip_p_b(P, S, 0); out$p0_guard <- r0$p; out$why0_guard <- r0$why   # test at 0 with the D-57 rules
  rt <- flip_p_b(P, S, theta); out$p_theta <- rt$p; out$theta_why <- rt$why
  Xf <- cbind(P$X0, P$x); k <- ncol(Xf)
  f <- tryCatch(suppressWarnings(glm.fit(Xf, P$y, family = binomial(), control = GLMC)), error = function(e) NULL)
  if (is.null(f) || !isTRUE(f$converged) || !all(is.finite(f$coefficients))) { out$est_why <- "full fit failed"; return(out) }
  if (max(abs(f$coefficients)) >= CI_DIVERGE) { out$est <- unname(f$coefficients[k]); out$est_why <- "full fit diverged"; return(out) }
  bh <- unname(f$coefficients[k]); out$est <- bh
  if (!(abs(bh) < CI_BMAX)) { out$est_why <- "|estimate| >= 10"; return(out) }
  w <- f$fitted.values * (1 - f$fitted.values)
  s <- tryCatch(sqrt(solve(crossprod(Xf, Xf * w))[k, k]), error = function(e) NA_real_, warning = function(w) NA_real_)
  out$se_mod <- s
  rc <- flip_p_b(P, S, bh); out$p_est <- rc$p; out$su_est <- rc$su
  L <- ci_side(P, S, bh, s, -1); U <- ci_side(P, S, bh, s, +1)
  out$lo <- L$end; out$hi <- U$end; out$rej_lo <- L$rej; out$rej_hi <- U$rej
  out$n_eval <- 3L + L$ne + U$ne; out$n_na <- L$nna + U$nna
  if (is.finite(L$end)) { out$p_lo <- flip_p_b(P, S, L$end)$p; out$p_lo_out <- flip_p_b(P, S, L$end - CI_EPS_BEYOND)$p }
  if (is.finite(U$end)) { out$p_hi <- flip_p_b(P, S, U$end)$p; out$p_hi_out <- flip_p_b(P, S, U$end + CI_EPS_BEYOND)$p }
  out
}
