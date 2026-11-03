# =============================================================================
# verify_plim.R -- independent test of the second review's central claim:
# that the GLMM's "coverage degrades as N grows" is largely an artefact of
# scoring `glmer` against 0.232, a value the estimator does not converge to.
#
# Test, using only stored fits (no refitting):
#  1. What does `glmer` actually converge to? Take the mean estimate at the
#     LARGEST N, where Monte Carlo noise is smallest -- that is the empirical
#     probability limit, and it needs no theory to compute.
#  2. Recompute coverage against BOTH targets and look at the trend in N. If the
#     review is right, the trend reverses; if it is wrong, it merely shifts.
#  3. Check the claim is specific: it should NOT hold at T = 2, and should be
#     weaker for GLMM-AR1.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr)})
RES <- path.expand("~/MMM/results_full")

load_cells <- function(pref, model, T, Ns) {
  pat <- sprintf("^%s_m%d_N(%s)_T%d_mmlb_c", pref, model, paste(Ns, collapse="|"), T)
  fs <- list.files(RES, pattern = pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  bind_rows(lapply(fs, readRDS)) %>%
    filter(term == "time:trt", method == "glmer", converged,
           is.finite(estimate), is.finite(se), abs(estimate) < 10, se > 0, se < 10)
}

NS <- c(25, 45, 95, 145)
for (spec in list(c("GenBinaryY", 1, 10), c("GenBinaryY", 1, 2), c("GLMMAR1", 1, 10))) {
  pref <- spec[1]; m <- as.integer(spec[2]); T <- as.integer(spec[3])
  d <- load_cells(pref, m, T, NS)
  if (is.null(d) || !nrow(d)) { cat(sprintf("\n%s M%d T=%d: no data\n", pref, m, T)); next }

  # empirical probability limit: the mean estimate at the largest available N
  big <- d %>% filter(N == max(N))
  plim <- mean(big$estimate)
  cat(sprintf("\n===== %s, model %d, T = %d =====\n", pref, m, T))
  cat(sprintf("mean estimate at N = %d: %.4f  (n = %d, MCSE %.4f)\n",
              max(d$N), plim, nrow(big), sd(big$estimate)/sqrt(nrow(big))))
  cat(sprintf("the analysis scores against 0.232 -> offset %.4f\n", plim - 0.232))

  z <- qnorm(0.975)
  tab <- d %>% group_by(N) %>%
    summarise(n = n(),
              cov_232  = 100*mean(abs(estimate - 0.232) <= z*se),
              cov_plim = 100*mean(abs(estimate - plim ) <= z*se),
              mcse     = 100*sqrt(0.9*0.1/n()), .groups = "drop")
  print(as.data.frame(tab %>% mutate(across(where(is.numeric), ~round(.x, 2)))),
        row.names = FALSE)
  tr232  <- tab$cov_232 [nrow(tab)] - tab$cov_232 [1]
  trplim <- tab$cov_plim[nrow(tab)] - tab$cov_plim[1]
  cat(sprintf("trend over N: vs 0.232 %+0.1f pp | vs the limit %+0.1f pp  -> %s\n",
              tr232, trplim,
              if (tr232 < -1 && trplim > 0) "REVERSES (review confirmed)"
              else if (tr232 < -1 && trplim > tr232 + 2) "shrinks but does not reverse"
              else "no artefact of this kind"))
}
