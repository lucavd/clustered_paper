# s7c_posthoc.R [workers] -- S7c, analyses added AFTER the pilot results (declared post hoc, exploratory).
# Why: B and C (residual time-trend r_i = sum_j time_ij e_ij) were far below 5% under H0 and had no power under H1.
# Diagnosis computed here: the per-subject r_i carries the subject's level residual (time-weighted), which cancels in the
# score S1 = sum_{trt=1} r_i because the reduced model forces sum_{trt=1} e_ij = 0, but inflates the between-subject
# variance used by the t-test. Variants:
#   B2 robust_t   robust score test: sum_i u_i / sqrt(sum_i u_i^2), u_i effective contributions, t with N - 2 df
#   B3 within_t   Student t (N - 2) on the within-subject residual trend w_i = sum_j (time_ij - mean_i(time)) e_ij
#   C2 within_perm permutation of arms among w_i, Welch statistic (as C)
#   A0 flip_relaxed  A without the |coefficient| < 10 filter on the reduced fit (converged, finite fitted values only)
# Diagnostics per replicate: m1 - m0 of r_i and its Student-t variance, for the variance-ratio explanation of B.
source(path.expand("~/MMM/R/s7c/s7c_common.R"))
a <- commandArgs(TRUE); W <- if (length(a)) as.integer(a[1]) else 20L
null_fit_relaxed <- function(P) {
  X0 <- cbind(1, P$time, P$trt)
  f <- tryCatch(suppressWarnings(glm.fit(X0, P$y, family = binomial(), control = glm.control(epsilon = 1e-12, maxit = 100))), error = function(e) NULL)
  if (is.null(f)) return(list(ok = FALSE))
  list(ok = isTRUE(f$converged) && all(is.finite(f$fitted.values)) && all(is.finite(f$coefficients)), mu = f$fitted.values, beta = f$coefficients, X0 = X0)
}
post_one <- function(env, N, ntime, label, i) {
  deg <- design_degenerate(1L, N, ntime, i); set_rep_seed(1L, i)
  d <- tryCatch(as.data.frame(gen_data(env, "M1", 1L, N, ntime, i, label)), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  P <- s7c_prep(d); g <- P$trt_s
  row <- data.table(sim_id = i, strict_ok = NA, relaxed_ok = FALSE, p_flip_relaxed = NA_real_, p_robust_t = NA_real_, p_within_t = NA_real_,
                    p_within_perm = NA_real_, dm_r = NA_real_, v_r = NA_real_, dm_w = NA_real_, v_w = NA_real_, S1 = NA_real_)
  F <- null_fit_relaxed(P); row$strict_ok <- isTRUE(F$ok) && max(abs(F$beta)) < COEF_MAX
  if (!isTRUE(F$ok)) return(row)
  S <- tryCatch(score_parts(P, F), error = function(e) NULL)
  if (is.null(S)) return(row)              # information matrix singular: score not computable
  row$relaxed_ok <- TRUE; e <- P$y - F$mu
  tbar <- ave(P$time, P$id); w <- as.numeric(rowsum((P$time - tbar) * e, P$id, reorder = TRUE)[, 1])
  set.seed(RESAMPLE_SEED(i)); row$p_flip_relaxed <- test_flip(S$u_eff)$p
  if (sum(S$u_eff^2) > 0) { tr <- sum(S$u_eff) / sqrt(sum(S$u_eff^2)); row$p_robust_t <- 2 * pt(-abs(tr), P$N - 2) }
  if (sum(g) >= 2 && sum(1 - g) >= 2) {
    st <- student_t(w, g); if (is.finite(st$t)) row$p_within_t <- 2 * pt(-abs(st$t), st$df)
    set.seed(RESAMPLE_SEED(i)); row$p_within_perm <- test_subj_perm(w, g)$p
    r <- S$r; n1 <- sum(g); n0 <- P$N - n1
    sp <- function(x) ((n1 - 1) * var(x[g == 1]) + (n0 - 1) * var(x[g == 0])) / (P$N - 2) * (1 / n1 + 1 / n0)
    row$dm_r <- mean(r[g == 1]) - mean(r[g == 0]); row$v_r <- sp(r)
    row$dm_w <- mean(w[g == 1]) - mean(w[g == 0]); row$v_w <- sp(w); row$S1 <- sum(r[g == 1])
  }
  row
}
tk <- pilot_tasks(); tk[, pfile := sub("/pilot/", "/posthoc/", file)]; dir.create(file.path(OUT_S7C, "posthoc"), showWarnings = FALSE)
todo <- tk[!file.exists(pfile)]
res <- parallel::mclapply(seq_len(nrow(todo)), function(k) tryCatch({
  t <- todo[k]; env <- load_env(t$dgp, RUNDIR_V4); ids <- (t$chunk - 1L) * 100L + seq_len(100L)
  R <- rbindlist(lapply(ids, function(i) post_one(env, t$N, t$ntime, t$label, i)), fill = TRUE)
  R[, `:=`(dgp = t$dgp, N = t$N, ntime = t$ntime, label = t$label)]; saveRDS(R, t$pfile); nrow(R) },
  error = function(e) paste("ERROR", conditionMessage(e))), mc.cores = W, mc.preschedule = FALSE)
cat("errors", sum(vapply(res, is.character, TRUE)), "\n")
