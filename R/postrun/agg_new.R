# =============================================================================
# agg_new.R -- aggregation for the redesigned scenarios, incorporating the
# three corrections the second review asked for. Reads stored chunks only.
#
# It runs on PARTIAL output by design: cells are complete or absent, never
# half-written, so aggregating mid-run exercises the whole path early instead
# of discovering a defect after eight hours.
#
# What is new relative to the previous aggregation:
#  (1) per-method exclusion filter. The likelihood-ratio p-value does not come
#      from the Wald standard error, so filtering those rows on `se` discards
#      valid results -- it previously threw away ~44,000 and understated that
#      method's level.
#  (2) coverage against BOTH targets: the generating parameter (conventional)
#      and the estimator's own empirical limit. The mixed model does not
#      converge to the value it was scored against, and scoring against a
#      non-estimand drives coverage to zero as N grows, by construction.
#  (3) both denominators: each method's own usable set, and the intersection
#      where every method returned a result. Cells with an EMPTY intersection
#      are counted and reported rather than averaged away -- under the old
#      design they were 59 of 480 and drove half of a manuscript claim.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr)})
RES <- path.expand("~/MMM/results_full")
OUT <- path.expand("~/MMM/tables")
dir.create(OUT, showWarnings = FALSE)

DGPS <- c(GenBinaryY = "GenBinaryY", GLMMnative = "GLMM-native", GLMMAR1 = "GLMM-AR1")
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink")
LRT  <- "glmer_lrt"
NMETH <- length(c(WALD, LRT))
z <- qnorm(0.975)

# Marginal truth per scenario. M4 deliberately uses a small effect; the others
# share a common target so the three mechanisms compare at equal effect size.
MARG <- c("1" = 0.200, "2" = 0.200, "3" = 0.200, "4" = 0.050, "5" = 0.200)

keep_rows <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method == LRT, converged, is.finite(p))     # no se filter here
)

load_complete <- function(pref, models) {
  fs <- list.files(RES, pattern = paste0("^", pref, "_m[", paste(models, collapse=""), "]_"),
                   full.names = TRUE)
  if (!length(fs)) return(NULL)
  key <- sub("_c[0-9]+[.]rds$", "", basename(fs))
  full <- names(which(table(key) == 10))          # only cells with all 10 chunks
  fs <- fs[key %in% full]
  if (!length(fs)) return(NULL)
  d <- bind_rows(lapply(fs, readRDS))
  d$model <- as.integer(d$model)
  d %>% filter(term == "time:trt")
}

LV <- list(); CV <- list(); PL <- list(); EM <- list(); k <- 0L
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  d <- load_complete(pref, 3:5)
  if (is.null(d)) { cat(sprintf("%-12s no complete cells yet\n", lab)); next }
  H0 <- keep_rows(d %>% filter(label == "manual"))
  H1 <- keep_rows(d %>% filter(label == "mmlb"))
  if (!nrow(H0) || !nrow(H1)) next

  # the mixed model's empirical limit, per cell shape, at the largest N present
  plim <- H1 %>% filter(method == "glmer") %>% group_by(model, ntime) %>%
    filter(N == max(N)) %>%
    summarise(plim = mean(estimate), n_plim = n(),
              mcse_plim = sd(estimate)/sqrt(n()), .groups = "drop") %>%
    mutate(dgp = lab)

  common <- function(df) df %>% distinct(model, N, ntime, sim_id, method) %>%
    count(model, N, ntime, sim_id, name = "nm") %>% filter(nm == NMETH) %>%
    select(model, N, ntime, sim_id)
  c0 <- common(H0); c1 <- common(H1)

  lvl <- function(df, s) df %>% group_by(method, model, N, ntime) %>%
    summarise(n = n(), k_rej = sum(p < 0.05), .groups = "drop") %>% mutate(set = s)

  tt <- expand.grid(model = sort(unique(H1$model)), ntime = sort(unique(H1$ntime)))
  tt$t_marg <- as.numeric(MARG[as.character(tt$model)])

  cov1 <- function(df, s) df %>% filter(method %in% WALD) %>%
    left_join(plim %>% select(model, ntime, plim), by = c("model","ntime")) %>%
    left_join(tt, by = c("model","ntime")) %>%
    mutate(t_gen  = t_marg,
           t_plim = ifelse(method == "glmer" & is.finite(plim), plim, t_marg)) %>%
    filter(is.finite(t_gen), is.finite(t_plim)) %>%
    mutate(cg = abs(estimate - t_gen) <= z*se,
           cp = abs(estimate - t_plim) <= z*se) %>%
    group_by(method, model, N, ntime) %>%
    summarise(n = n(), cov_gen = 100*mean(cg), cov_plim = 100*mean(cp),
              mcse = 100*sqrt(mean(cg)*(1-mean(cg))/n()), .groups = "drop") %>%
    mutate(set = s)

  k <- k + 1L
  LV[[k]] <- bind_rows(lvl(H0, "own"),
                       lvl(H0 %>% inner_join(c0, by = c("model","N","ntime","sim_id")),
                           "intersection")) %>% mutate(dgp = lab)
  CV[[k]] <- bind_rows(cov1(H1, "own"),
                       cov1(H1 %>% inner_join(c1, by = c("model","N","ntime","sim_id")),
                            "intersection")) %>% mutate(dgp = lab)
  PL[[k]] <- plim
  EM[[k]] <- H0 %>% distinct(model, N, ntime) %>%
    anti_join(c0 %>% distinct(model, N, ntime), by = c("model","N","ntime")) %>%
    mutate(dgp = lab)
  rm(d, H0, H1); gc(verbose = FALSE)
}
if (!length(LV)) { cat("\nnot enough complete cells yet\n"); quit(status = 0) }
L <- bind_rows(LV); C <- bind_rows(CV); P <- bind_rows(PL); E <- bind_rows(EM)
write.csv(L, file.path(OUT, "N_level_m35.csv"),    row.names = FALSE)
write.csv(C, file.path(OUT, "N_coverage_m35.csv"), row.names = FALSE)
write.csv(P, file.path(OUT, "N_glmer_plim.csv"),   row.names = FALSE)
write.csv(E, file.path(OUT, "N_empty_cells.csv"),  row.names = FALSE)

cat(sprintf("\ncells contributing: %d | replicate rows: %s\n",
            n_distinct(L$model, L$N, L$ntime, L$dgp), format(sum(L$n), big.mark = ",")))

cat("\n===== does the mixed model now FIT where it could not before? =====\n")
est <- L %>% filter(set == "own", method %in% c("glmer","glmer_lrt","mmm")) %>%
  group_by(dgp, model, ntime) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = ntime, values_from = n, names_prefix = "T")
print(as.data.frame(est), row.names = FALSE)

cat("\n===== TYPE I ERROR, both denominators (nominal 5%) =====\n")
p1 <- L %>% group_by(dgp, method, set) %>%
  summarise(rate = 100*sum(k_rej)/sum(n), n = sum(n), .groups = "drop") %>%
  mutate(rate = round(rate, 2)) %>%
  pivot_wider(names_from = set, values_from = c(rate, n))
print(as.data.frame(p1 %>% arrange(dgp, rate_own)), row.names = FALSE)

cat("\n===== COVERAGE, both targets (nominal 95%) =====\n")
p2 <- C %>% filter(set == "own") %>% group_by(dgp, method) %>%
  summarise(gen = round(sum(cov_gen*n)/sum(n), 2),
            plim = round(sum(cov_plim*n)/sum(n), 2), n = sum(n), .groups = "drop")
print(as.data.frame(p2 %>% arrange(dgp, desc(gen))), row.names = FALSE)

cat("\n===== the mixed model's target vs where it converges =====\n")
print(as.data.frame(P %>% mutate(across(c(plim, mcse_plim), ~round(.x, 4))) %>%
                      arrange(dgp, model, ntime)), row.names = FALSE)

cat("\n===== cells with an EMPTY intersection =====\n")
if (nrow(E)) print(as.data.frame(E %>% count(dgp, model, ntime, name = "cells")),
                   row.names = FALSE) else cat("none\n")

cat("\nwrote N_level_m35.csv, N_coverage_m35.csv, N_glmer_plim.csv, N_empty_cells.csv\n")
