# =============================================================================
# s7d_gate_formulas.R -- gate GE4 of S7d (independent check of the formulas), run before the production run.
# (a) fewclust 0.2.0 gee_vcov md / kc / fg against a hand-coded implementation written from the articles, on the same
#     geeglm fit, for the three working correlations: Mancl-DeRouen (I - H_ii)^{-1} e_i; Kauermann-Carroll with the
#     principal inverse square root of the non-symmetric (I - H_ii), computed as Q diag(l^{-1/2}) Q^{-1}; Fay-Graubard
#     with diag{(1 - min(0.75, Q_ii,jj))^{-1/2}}, Q_ii = D_i' V_i^{-1} D_i B^{-1}. H_ii = D_i B^{-1} D_i' V_i^{-1}.
#     Pass: relative difference < 1e-8 wherever both are computed; datasets computed by only one of the two <= 1% per correlation.
# (b) sign-flip score in the scenarios with a covariate (M2 age, M5 dose): sum of the effective contributions equals the
#     numDeriv gradient of the full-model log-likelihood at (reduced fit, time:trt = 0); relative difference < 1e-6.
# Datasets: every scenario x mechanism (22), N in {10, 25}, T in {3, 10}, i = 1..5, H0 and H1 (880 datasets).
# Also recorded (not gate): geesmv's KC formula (eigenvectors transposed instead of inverted) on the same fit.
# =============================================================================
source(path.expand("~/MMM/R/s7d/s7d_common.R"))
suppressMessages(library(numDeriv))
hand <- function(fit) {
  X <- model.matrix(fit); mu <- fit$fitted.values; a <- mu * (1 - mu); y <- fit$y
  al <- unname(fit$geese$alpha); cs <- fit$corstr; ids <- split(seq_along(y), factor(fit$id, levels = unique(fit$id)))
  Rm <- function(n) switch(cs, independence = diag(n), exchangeable = { m <- matrix(al, n, n); diag(m) <- 1; m }, ar1 = al^abs(outer(1:n, 1:n, "-")))
  cl <- lapply(ids, function(ix) { n <- length(ix); Di <- X[ix, , drop = FALSE] * a[ix]; s <- sqrt(a[ix]); Vi <- (s %o% s) * Rm(n)
    list(D = Di, Vinv = solve(Vi), e = y[ix] - mu[ix]) })
  B <- Reduce(`+`, lapply(cl, function(z) t(z$D) %*% z$Vinv %*% z$D)); Bi <- solve(B)
  k <- match("time:trt", colnames(X))
  meat <- function(f) { M <- Reduce(`+`, lapply(cl, function(z) { u <- f(z); u %*% t(u) })); sqrt((Bi %*% M %*% Bi)[k, k]) }
  IH <- function(z) diag(nrow(z$D)) - z$D %*% Bi %*% t(z$D) %*% z$Vinv
  pr <- function(A) { E <- eigen(A); if (any(Mod(E$values) <= 1e-8) || any(Re(E$values) <= 0)) stop("leverage"); Re(E$vectors %*% diag(1 / sqrt(E$values), length(E$values)) %*% solve(E$vectors)) }
  gs <- function(A) { E <- eigen(A); d <- Re(E$values); d <- (d + abs(d)) / 2; d2 <- 1 / sqrt(d); d2[d == 0] <- 0; Re(E$vectors %*% diag(d2, length(d)) %*% t(E$vectors)) }
  f1 <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  c(md = f1(meat(function(z) t(z$D) %*% z$Vinv %*% solve(IH(z), z$e))),
    kc = f1(meat(function(z) t(z$D) %*% z$Vinv %*% pr(IH(z)) %*% z$e)),
    fg = f1(meat(function(z) { Q <- t(z$D) %*% z$Vinv %*% z$D %*% Bi; (t(z$D) %*% z$Vinv %*% z$e) / sqrt(1 - pmin(0.75, diag(Q))) })),
    kc_geesmv_formula = f1(meat(function(z) t(z$D) %*% z$Vinv %*% gs(IH(z)) %*% z$e)))
}
TK <- unique(s7d_tasks()[, .(dgp, scen, model)])
G <- TK[, CJ(N = c(10L, 25L), ntime = c(3L, 10L), label = c("manual", "mmlb"), i = 1:5), by = .(dgp, scen, model)]
one <- function(g) {
  env <- load_env(g$dgp, RUNDIR_V4); gen <- s7d_gen(env, g)
  set_rep_seed(g$model, g$i); d <- tryCatch(as.data.frame(gen(g$N, g$ntime, g$i, g$label)), error = function(e) NULL)
  if (is.null(d)) return(NULL); d <- d[order(d$id, d$time), ]; form <- FORM_OF(g$scen); out <- list()
  for (cn in names(CORS)) {
    fit <- tryCatch(geepack::geeglm(form, id = id, data = d, family = binomial, corstr = CORS[[cn]], scale.fix = TRUE), error = function(e) NULL)
    if (is.null(fit) || !all(is.finite(coef(fit))) || max(abs(coef(fit))) >= 10) next
    k <- match("time:trt", names(coef(fit))); h <- tryCatch(hand(fit), error = function(e) c(md = NA, kc = NA, fg = NA, kc_geesmv_formula = NA))
    fc <- sapply(c("md", "kc", "fg"), function(ty) tryCatch(sqrt(fewclust::gee_vcov(fit, ty)[k, k]), error = function(e) NA_real_))
    out[[cn]] <- data.table(g, cor = cn, h_md = h[["md"]], h_kc = h[["kc"]], h_fg = h[["fg"]], h_kc_geesmv = h[["kc_geesmv_formula"]],
                            fc_md = fc[["md"]], fc_kc = fc[["kc"]], fc_fg = fc[["fg"]])
  }
  R <- rbindlist(out)
  sc <- NULL
  if (g$scen %in% c("M2", "M5")) {
    id <- as.integer(factor(d$id)); y <- d$binY; x <- d$time * d$trt; X0 <- s7d_reduced_X(d, g$scen)
    f <- suppressWarnings(glm.fit(X0, y, family = binomial(), control = glm.control(epsilon = 1e-12, maxit = 100)))
    if (isTRUE(f$converged) && all(is.finite(f$coefficients))) {
      mu <- f$fitted.values; e <- y - mu; w <- mu * (1 - mu)
      u <- as.numeric(rowsum(x * e, id, reorder = TRUE)[, 1] - rowsum(X0 * e, id, reorder = TRUE) %*% solve(crossprod(X0, X0 * w), colSums(X0 * (w * x))))
      ll <- function(b) { eta <- cbind(X0, x) %*% b; sum(y * eta - log1p(exp(eta))) }
      nd <- numDeriv::grad(ll, c(f$coefficients, 0))[ncol(X0) + 1]
      sc <- data.table(g, sum_u = sum(u), numderiv = nd)
    }
  }
  list(R = R, sc = sc)
}
res <- parallel::mclapply(seq_len(nrow(G)), function(r) tryCatch(one(G[r]), error = function(e) NULL), mc.cores = 20L, mc.preschedule = FALSE)
R <- rbindlist(lapply(res, `[[`, "R"), fill = TRUE); SC <- rbindlist(lapply(res, `[[`, "sc"), fill = TRUE)
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)
chk <- R[, .(n = .N, both_md = sum(is.finite(h_md) & is.finite(fc_md)), one_md = sum(xor(is.finite(h_md), is.finite(fc_md))), max_md = max(rel(fc_md, h_md), na.rm = TRUE),
             both_kc = sum(is.finite(h_kc) & is.finite(fc_kc)), one_kc = sum(xor(is.finite(h_kc), is.finite(fc_kc))), max_kc = max(rel(fc_kc, h_kc), na.rm = TRUE),
             both_fg = sum(is.finite(h_fg) & is.finite(fc_fg)), one_fg = sum(xor(is.finite(h_fg), is.finite(fc_fg))), max_fg = max(rel(fc_fg, h_fg), na.rm = TRUE)), by = cor]
print(chk)
sc_max <- if (nrow(SC)) max(rel(SC$sum_u, SC$numderiv)) else NA_real_
cat("flip score vs numDeriv: datasets", nrow(SC), "max rel diff", format(sc_max, digits = 3), "\n")
kcg <- R[is.finite(h_kc_geesmv) & is.finite(h_kc), .(n = .N, med_ratio = median(h_kc_geesmv / h_kc), q05 = quantile(h_kc_geesmv / h_kc, 0.05), q95 = quantile(h_kc_geesmv / h_kc, 0.95),
                                                     share_off_1pct = mean(abs(h_kc_geesmv / h_kc - 1) > 0.01)), by = cor]
print(kcg)
pass_a <- all(c(chk$one_md, chk$one_kc, chk$one_fg) <= 0.01 * rep(chk$n, 3)) && all(c(chk$max_md, chk$max_kc, chk$max_fg) < 1e-8)
pass_b <- nrow(SC) >= 50 && sc_max < 1e-6
dir.create(file.path(OUT_S7D, "gates"), showWarnings = FALSE)
fwrite(R, file.path(OUT_S7D, "gates/GE4_formulas.csv")); fwrite(SC, file.path(OUT_S7D, "gates/GE4_flip_score.csv"))
fwrite(chk, file.path(REPO, "tables/V4_S7d_GE4.csv")); fwrite(kcg, file.path(REPO, "tables/V4_S7d_kc_geesmv_formula.csv"))
cat("GE4", if (pass_a && pass_b) "PASS" else "FAIL", "(a)", pass_a, "(b)", pass_b, "\n")
