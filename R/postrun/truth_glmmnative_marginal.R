# =============================================================================
# truth_ar1.R -- D1/D2: recover the TRUE parameter values each mechanism used,
# instead of the hard-coded constants the metrics module assumes.
#
# D1: the AR(1) generator solves for its conditional coefficients per replicate
#     against total latent variance D + sigma_ar^2, but nothing ever read the
#     result back -- performance_metrics.R hard-codes cond = 0.232 for all three
#     mechanisms.
# D2: the GLMM-native arm sets beta_c = 0.2 * 1.16006 (the Zeger approximation).
#     The marginal interaction that ACTUALLY implies is not 0.2, and it drifts
#     with T.
#
# No model is fitted here. The design is a deterministic function of the seed,
# so we regenerate the design only, recover each mechanism's own truth, and
# write a per-cell truth table for the aggregation to join on.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(fastGHQuad)})
MSM <- path.expand("~/MMM/MMM_sim_manual")
sys.source(file.path(MSM, "gen_ar1_data.R"), envir = globalenv())

N_LEVELS <- c(10, 15, 20, seq(25, 145, 10)); T_LEVELS <- c(2, 3, 4, 6, 8, 10)
MODELS <- 1:5

# ---- exact marginalisation of a conditional logistic over a normal intercept ----
gh <- fastGHQuad::gaussHermiteData(80)
marg_p <- function(eta_c, D) sum(gh$w * plogis(eta_c + sqrt(2*D)*gh$x)) / sqrt(pi)

# population independence-GEE limit: fit a marginal logistic to the EXACT
# marginal probabilities implied by the conditional model over the design
implied_marginal <- function(bc_int, T, b0, btime, btrt, D) {
  g <- expand.grid(time = 0:(T-1), trt = 0:1)
  X <- model.matrix(~ time*trt, g)
  eta_c <- b0 + btime*g$time + btrt*g$trt + bc_int*g$time*g$trt
  p <- vapply(eta_c, marg_p, numeric(1), D = D)
  beta <- rep(0, ncol(X))
  for (it in 1:200) {
    mu <- plogis(X %*% beta); W <- as.vector(mu*(1-mu))
    z  <- X %*% beta + (p - mu)/pmax(W, 1e-12)
    nb <- solve(crossprod(X, X*W), crossprod(X, W*z))
    if (max(abs(nb - beta)) < 1e-12) { beta <- nb; break }; beta <- nb
  }
  as.vector(beta)[4]
}

cat("===== D2: GLMM-native -- what marginal effect does beta_c = 0.232 imply? =====\n")
d2 <- do.call(rbind, lapply(T_LEVELS, function(T) {
  data.frame(ntime = T,
             m_main = implied_marginal(0.232, T, -2.32, 0.116, 0.232, 1),
             m_m4   = implied_marginal(0.058, T, -2.32, 0.232, 0.232, 1))
}))
d2$dev_main <- d2$m_main - 0.20; d2$dev_m4 <- d2$m_m4 - 0.05
print(as.data.frame(d2 %>% mutate(across(where(is.numeric), ~round(.x, 5)))), row.names = FALSE)
cat("\n=> the analysis scored this arm against a flat 0.200 / 0.050.\n")

cat("\n===== D1: AR(1) -- what conditional value did the generator actually use? =====\n")
set.seed(20260808)
grid <- expand.grid(model = MODELS, ntime = T_LEVELS, N = c(25, 75, 145))
out <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  m <- grid$model[k]; T <- grid$ntime[k]; N <- grid$N[k]
  vals <- numeric(0)
  for (rep in 1:12) {
    bc <- tryCatch(calibrate_beta_c_ar1(N = N, ntime = T, model = m), error = function(e) NULL)
    if (!is.null(bc)) vals <- c(vals, bc[length(bc)])
  }
  if (length(vals))
    out[[k]] <- data.frame(model = m, ntime = T, N = N,
                           cond_used = mean(vals), sd = sd(vals), n_ok = length(vals))
}
A <- bind_rows(out)
if (nrow(A)) {
  cat("mean conditional interaction actually used, by model and T:\n")
  print(as.data.frame(A %>% group_by(model, ntime) %>%
    summarise(cond_used = round(mean(cond_used), 4), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = ntime, values_from = cond_used,
                       names_prefix = "T=")), row.names = FALSE)
  cat("\nthe metrics module hard-codes 0.232 (model 4: 0.058) for every cell.\n")
  write.csv(A, path.expand("~/MMM/tables/F_truth_ar1.csv"), row.names = FALSE)
} else cat("calibrate_beta_c_ar1() not callable with this signature -- see log\n")
write.csv(d2, path.expand("~/MMM/tables/F_truth_glmmnative.csv"), row.names = FALSE)
cat("\nwrote F_truth_ar1.csv, F_truth_glmmnative.csv\n")
