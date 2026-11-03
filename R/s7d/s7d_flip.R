# =============================================================================
# s7d_flip.R -- sign-flip score test (D-46, D-47) generalised to the mean model of each scenario.
# Reduced model = the scenario's GEE mean model without time:trt (M2: + age, M5: + dose; else time + trt),
# independence working correlation, glm.fit(epsilon = 1e-12, maxit = 100). For M1/M1c/M3/M4/M6/M7 this is
# exactly R/s7c/s7c_tests.R (X0 = (1, time, trt)). Effective score contributions per subject (Hemerik,
# Goeman, Finos 2020); test_flip() from R/s7c/s7c_tests.R (all 2^N signs if N <= 12, else 1,999 + identity).
# Two versions: "flip" (primary, D-47: converged, finite coefficients and fitted values, information invertible)
# and "flip_strict" (sensitivity: in addition max |coefficient| < 10, the pre-registered S7c filter).
# =============================================================================
s7d_reduced_X <- function(d, scen) {
  base <- cbind(1, as.numeric(d$time), as.numeric(d$trt))
  if (scen == "M2") base <- cbind(base, as.numeric(d$age))
  if (scen == "M5") base <- cbind(base, as.numeric(d$dose))
  base
}
s7d_flip <- function(d, scen) {   # d sorted by id, time; RNG state set by the caller
  id <- as.integer(factor(d$id)); y <- as.numeric(d$binY); x <- as.numeric(d$time) * as.numeric(d$trt)
  X0 <- s7d_reduced_X(d, scen)
  f <- tryCatch(suppressWarnings(glm.fit(X0, y, family = binomial(), control = glm.control(epsilon = 1e-12, maxit = 100))),
                error = function(e) NULL)
  na <- data.frame(test = c("flip", "flip_strict"), p = NA_real_, why = NA_character_, max_abs_coef = NA_real_, n_ref = NA_integer_)
  if (is.null(f)) { na$why <- "reduced fit error"; return(na) }
  b <- f$coefficients; mab <- if (all(is.finite(b))) max(abs(b)) else NA_real_; na$max_abs_coef <- mab
  if (!(isTRUE(f$converged) && all(is.finite(b)) && all(is.finite(f$fitted.values)))) {
    na$why <- if (!isTRUE(f$converged)) "reduced fit not converged" else "reduced fit non-finite"; return(na) }
  mu <- f$fitted.values; e <- y - mu; w <- mu * (1 - mu)
  I00 <- crossprod(X0, X0 * w); Ix0 <- colSums(X0 * (w * x))
  proj <- tryCatch(solve(I00, Ix0), error = function(e) NULL)
  if (is.null(proj)) { na$why <- "information singular"; return(na) }
  u <- as.numeric(rowsum(x * e, id, reorder = TRUE)[, 1] - rowsum(X0 * e, id, reorder = TRUE) %*% proj)
  a <- test_flip(u)
  out <- na; out$p[1] <- a$p; out$why[1] <- a$why; out$n_ref[1] <- a$n_ref
  if (mab < COEF_MAX) { out$p[2] <- a$p; out$why[2] <- a$why; out$n_ref[2] <- a$n_ref } else out$why[2] <- "reduced fit coefficient >= 10"
  out
}
