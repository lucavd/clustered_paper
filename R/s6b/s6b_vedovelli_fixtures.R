# =============================================================================
# s6b_vedovelli_fixtures.R -- test fixtures for a future implementation of the "Vedovelli procedure" for MMM
# (convergence check of D-32, D-35) in the fewclust package. Picks study datasets of M3 (GenBinaryY, N = 25,
# T = 8, H0) covering every branch of the rule, regenerates them with their seed (as R/s6/s6_mmm_cc_run.R),
# refits A and B to confirm the stored values, and writes the data and the expected results.
# Output: ~/MMM-package/fewclust/dev/vedovelli_fixtures/{case_*.csv, expected.csv, README.md}
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
suppressMessages({ library(MMLB); library(dplyr) })
RD <- file.path(REPO, "R/v4/rundir")
for (f in c(sprintf("generation_model%d.R", 1:5), "generate_data.R", "make_all_data.R")) source(file.path(RD, f))
OUT <- path.expand("~/MMM-package/fewclust/dev/vedovelli_fixtures"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CC <- readRDS(file.path(REPO, "results_v4/S6cc/cc_all.rds"))
X <- CC[dgp == "GenBinaryY" & model == 3 & label == "manual" & A_ok %in% TRUE & A_has_term %in% TRUE & N == 25 & ntime == 8]
X[, branch := fifelse(!bnd, "keep_converged", fifelse(!(A_code %in% c(1L, 2L)) & used == "B", "fallback_nonconvergence",
              fifelse(used == "B", "fallback_slope_at_boundary", "fallback_failed_keep_study_fit")))]
print(X[, .N, by = branch])
pick <- X[, head(.SD[order(sim_id)], 2), by = branch]
fitw <- function(d, lv) { x <- tryCatch(suppressWarnings(mm(binY ~ time * trt, lv.formula = lv, t.formula = ~1, data = d, id = id, q = 20)), error = function(e) NULL)
  if (is.null(x)) return(c(est = NA, se = NA, code = NA)); s <- summary(x)$mean.table; k <- which(rownames(s) == "time:trt")
  c(est = s[k, 1], se = s[k, 2], code = unname(x$control["convergence_code"])) }
same <- function(a, b) (is.finite(a) & is.finite(b) & abs(a - b) < 1e-8) | (!is.finite(a) & !is.finite(b))
res <- list()
for (r in seq_len(nrow(pick))) {
  k <- pick[r]
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(20260805L + 1000L * 3L + k$sim_id)
  d <- as.data.frame(make_all_data(k$N, k$ntime, k$sim_id, 3L, k$label))
  fn <- sprintf("case_%02d_%s.csv", r, k$branch); fwrite(d[, c("id", "time", "trt", "binY")], file.path(OUT, fn))
  A <- fitw(d, ~time); B <- if (k$bnd) fitw(d, ~1) else c(est = NA, se = NA, code = NA)
  res[[r]] <- data.table(file = fn, branch = k$branch, sim_id = k$sim_id, n_clusters = k$N, A_est = k$A_est, A_se = k$A_se, A_code = k$A_code,
    A_latent_sd_slope = k$A_slope_min, at_boundary = k$bnd, used = k$used, est = k$cc_est, se = k$cc_se,
    p_normal = k$cc_p, p_t = 2 * pt(-abs(k$cc_est / k$cc_se), df = k$N - 4L), df_t = k$N - 4L,
    refit_A_matches = same(A[["est"]], k$A_est) & same(A[["se"]], k$A_se), refit_B_matches = if (k$used == "B") same(B[["est"]], k$cc_est) & same(B[["se"]], k$cc_se) else NA)
}
E <- rbindlist(res); fwrite(E, file.path(OUT, "expected.csv")); print(E[, .(file, used, refit_A_matches, refit_B_matches)])
stopifnot(all(E$refit_A_matches), all(E$refit_B_matches %in% c(TRUE, NA)))
