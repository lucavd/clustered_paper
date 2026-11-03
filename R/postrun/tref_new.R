# =============================================================================
# tref_new.R -- recompute the reference-distribution comparison on the corrected
# design. Section 3.7 still carried numbers from the pre-correction grid.
#
# The Wald statistic is referred to a t distribution with between-within degrees
# of freedom (N - p) instead of the normal. This acts on the test statistic, not
# on the mean model, so the direction of the effect cannot change -- but the
# magnitudes must be recomputed, and which method benefits most is an empirical
# question that the corrected scenarios can answer differently.
#
# The likelihood-ratio test is excluded: its p-value comes from a chi-square on
# the deviance difference and no Wald statistic is involved.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/tables")

DGPS <- c(GenBinaryY = "GenBinaryY", GLMMnative = "GLMM-native", GLMMAR1 = "GLMM-AR1")
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink")
# free parameters in the mean model, per scenario (used for the df)
NP <- c("1" = 4, "2" = 5, "3" = 4, "4" = 4, "5" = 5)

out <- list(); k <- 0L
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs  <- list.files(RES, pattern = paste0("^", pref, "_"), full.names = TRUE)
  key <- sub("_c[0-9]+[.]rds$", "", basename(fs))
  fs  <- fs[key %in% names(which(table(key) == 10))]
  d <- bind_rows(lapply(fs, readRDS))
  d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt", label == "manual", method %in% WALD,
                    converged, is.finite(estimate), is.finite(se),
                    abs(estimate) < 10, se > 0, se < 10)
  d$df <- d$N - NP[as.character(d$model)]
  d <- d %>% filter(df > 0) %>%
    mutate(z = estimate / se,
           rej_norm = 2 * pnorm(-abs(z)) < 0.05,
           rej_t    = 2 * pt(-abs(z), df = df) < 0.05)
  k <- k + 1L
  out[[k]] <- d %>% group_by(method, N) %>%
    summarise(n = n(), norm = 100 * mean(rej_norm), tref = 100 * mean(rej_t),
              .groups = "drop") %>% mutate(dgp = lab)
  rm(d); gc(verbose = FALSE)
}
R <- bind_rows(out)
write.csv(R, file.path(OUT, "O_reference_dist.csv"), row.names = FALSE)

cat("===== pooled rejection rate, normal vs t reference (nominal 5%) =====\n")
p <- R %>% group_by(dgp, method) %>%
  summarise(norm = sum(norm * n) / sum(n), tref = sum(tref * n) / sum(n),
            n = sum(n), .groups = "drop") %>%
  mutate(shift = tref - norm) %>%
  mutate(across(c(norm, tref, shift), ~round(.x, 2)))
print(as.data.frame(p %>% arrange(dgp, tref)), row.names = FALSE)

cat("\n===== worst-case |deviation| across N, under each reference =====\n")
w <- R %>% group_by(dgp, method) %>%
  summarise(wc_norm = max(abs(norm - 5)), wc_t = max(abs(tref - 5)), .groups = "drop") %>%
  mutate(across(c(wc_norm, wc_t), ~round(.x, 2)))
print(as.data.frame(w %>% arrange(dgp, wc_t)), row.names = FALSE)

cat("\n===== who has the smallest worst case under each reference? =====\n")
best <- w %>% group_by(dgp) %>%
  summarise(best_normal = method[which.min(wc_norm)], val_n = min(wc_norm),
            best_t      = method[which.min(wc_t)],    val_t = min(wc_t), .groups = "drop")
print(as.data.frame(best), row.names = FALSE)

cat("\n===== the effect is largest at few clusters (N <= 25) =====\n")
s <- R %>% filter(N <= 25) %>% group_by(method) %>%
  summarise(norm = sum(norm * n) / sum(n), tref = sum(tref * n) / sum(n), .groups = "drop") %>%
  mutate(shift = round(tref - norm, 2), across(c(norm, tref), ~round(.x, 2))) %>%
  arrange(shift)
print(as.data.frame(s), row.names = FALSE)

cat("\nwrote O_reference_dist.csv\n")
