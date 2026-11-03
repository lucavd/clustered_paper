# s7c_gates_pre.R -- GD4 (formula cross-checks on 50 datasets) + flipscores cross-check (reported, not a gate)
source(path.expand("~/MMM/R/s7c/s7c_common.R"))
suppressMessages(library(numDeriv))
.libPaths(c(.libPaths(), S7C_LIB))
has_fs <- requireNamespace("flipscores", quietly = TRUE)
envs <- setNames(lapply(DGPS, load_env, rundir = RUNDIR_V4), DGPS)
grid <- CJ(dgp = DGPS, N = PILOT_N, ntime = PILOT_T, i = 1:4)[order(i)]
out <- list(); k <- 0
for (r in seq_len(nrow(grid))) {
  if (k >= 50) break
  g <- grid[r]; set_rep_seed(1L, g$i)
  d <- tryCatch(as.data.frame(gen_data(envs[[g$dgp]], "M1", 1L, g$N, g$ntime, g$i, "manual")), error = function(e) NULL)
  if (is.null(d)) next
  P <- s7c_prep(d); F <- null_fit(P); if (!F$ok) next
  S <- score_parts(P, F); x <- P$time * P$trt; eta0 <- as.numeric(F$X0 %*% F$beta)
  ll <- function(b3) { eta <- eta0 + b3 * x; sum(P$y * eta - log1p(exp(eta))) }
  gr <- numDeriv::grad(ll, 0)
  g1 <- P$trt_s; n1 <- sum(g1); n0 <- P$N - n1
  if (n1 < 2 || n0 < 2) next
  k <- k + 1
  st <- student_t(S$r, g1); tt <- t.test(S$r[g1 == 1], S$r[g1 == 0], var.equal = TRUE)$statistic
  wt <- t.test(S$r[g1 == 1], S$r[g1 == 0])$statistic; wm <- welch_mat(S$r, matrix(g1, 1))
  S1 <- sum(S$r[g1 == 1]); id_lhs <- S1; id_rhs <- n1 * n0 / P$N * (mean(S$r[g1 == 1]) - mean(S$r[g1 == 0]))
  fs_p <- NA_real_; my_p <- NA_real_; fs_msg <- NA_character_
  if (has_fs && P$N <= ENUM_FLIP_MAXN) {
    my_p <- test_flip(S$u_eff)$p
    fs <- tryCatch({ m <- flipscores::flipscores(binY ~ time + trt + time:trt, family = binomial, data = d,
                                                 score_type = "effective", id = d$id, n_flips = 2^P$N)
                     summary(m)$coefficients }, error = function(e) e)
    if (inherits(fs, "error")) fs_msg <- conditionMessage(fs) else { pc <- grep("^Pr", colnames(fs)); fs_p <- fs["time:trt", pc[length(pc)]] }
  }
  out[[k]] <- data.table(dgp = g$dgp, N = g$N, ntime = g$ntime, i = g$i,
    score_eff = sum(S$u_eff), score_raw = sum(S$u_raw), grad = gr,
    rel_eff = abs(sum(S$u_eff) - gr) / max(1, abs(gr)), rel_raw = abs(sum(S$u_raw) - gr) / max(1, abs(gr)),
    sum_r = sum(S$r), ident_err = abs(id_lhs - id_rhs), ident_exact_err = abs(id_lhs - n1 * sum(S$r) / P$N - id_rhs), null_conv = F$converged, t_err = abs(st$t - tt), welch_err = abs(wm - wt),
    flip_p = my_p, flipscores_p = fs_p, flipscores_msg = fs_msg)
}
G <- rbindlist(out); fwrite(G, file.path(OUT_S7C, "gd4_checks.csv"))
pass <- nrow(G) == 50 && all(G$rel_eff < 1e-6) && all(G$ident_err < 1e-10) && all(G$t_err < 1e-10) && all(G$welch_err < 1e-10)
cat("GD4 datasets", nrow(G), "max rel_eff", signif(max(G$rel_eff), 3), "max rel_raw", signif(max(G$rel_raw), 3),
    "max ident", signif(max(G$ident_err), 3), "max t", signif(max(G$t_err), 3), "max welch", signif(max(G$welch_err), 3), "PASS", pass, "\n")
cat("flipscores available", has_fs, "; compared", sum(is.finite(G$flipscores_p)), "; max |p diff|",
    if (any(is.finite(G$flipscores_p))) signif(max(abs(G$flip_p - G$flipscores_p), na.rm = TRUE), 3) else NA,
    "; errors", sum(!is.na(G$flipscores_msg)), if (any(!is.na(G$flipscores_msg))) G$flipscores_msg[!is.na(G$flipscores_msg)][1] else "", "\n")
