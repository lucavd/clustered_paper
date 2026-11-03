# =============================================================================
# reagg_all.R -- the three re-aggregations the second review asked for.
# Reads the stored chunks only; no simulation is re-run.
#
#  (1) Regenerate the bootstrap / by-scenario / paired tables WITH the
#      scenario-5 GEE refit swapped in. The committed F_* tables were built
#      before results_m5gee/ existed and still carry the defective rows.
#  (2) Score the mixed model's coverage against BOTH targets: the generating
#      parameter (conventional) and the estimator's own empirical limit. An
#      estimator scored against a value it does not converge to has coverage
#      going to zero as N grows, by construction, so the trend in N cannot
#      evidence "intervals too narrow".
#  (3) Add the intersection denominator to the coverage table, as was already
#      done for the level table, and report per-cell so the empty cells are
#      visible instead of being averaged away.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr)})
RES  <- path.expand("~/MMM/results_full")
M5   <- path.expand("~/MMM/results_m5gee")
OUT  <- path.expand("~/MMM/tables")
set.seed(20260808)

DGPS <- c(GenBinaryY = "GenBinaryY", GLMMnative = "GLMM-native", GLMMAR1 = "GLMM-AR1")
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink")
LRT  <- "glmer_lrt"
MARG <- c("1" = 0.200, "2" = 0.200, "3" = 0.200, "4" = 0.050, "5" = 0.200)

ar1 <- read.csv(file.path(OUT, "F_truth_ar1.csv")) %>%
  group_by(model, ntime) %>% summarise(cond = mean(cond_used), .groups = "drop")
gn  <- read.csv(file.path(OUT, "F_truth_glmmnative.csv"))

cond_gen <- function(dgp, model, ntime) {
  if (dgp == "GLMM-AR1") {
    v <- ar1$cond[ar1$model == model & ar1$ntime == ntime]
    if (length(v)) return(v[1])
  }
  if (model == 4) 0.058 else 0.232
}
marg_gen <- function(dgp, model, ntime) {
  if (dgp == "GLMM-native") {
    r <- gn[gn$ntime == ntime, ]
    if (nrow(r)) return(if (model == 4) r$marg_m4[1] else r$marg_main[1])
  }
  unname(MARG[as.character(model)])
}

keep_rows <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method == LRT, converged, is.finite(p))
)

z <- qnorm(0.975)
CW <- list(); CV <- list(); PL <- list(); SC <- list(); k <- 0L

for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs  <- list.files(RES, pattern = paste0("^", pref, "_"), full.names = TRUE)
  key <- sub("_c[0-9]+[.]rds$", "", basename(fs))
  d   <- bind_rows(lapply(fs[key %in% names(which(table(key) == 10))], readRDS))
  d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt")

  # (1) scenario-5 GEE refit
  m5f <- list.files(M5, pattern = paste0("^", pref, "_"), full.names = TRUE)
  if (length(m5f)) {
    m5 <- bind_rows(lapply(m5f, readRDS)) %>% filter(term == "time:trt")
    m5$model <- 5L
    d <- bind_rows(d %>% filter(!(model == 5L & method == "gee")),
                   m5[, intersect(names(m5), names(d)), drop = FALSE])
  }

  H0 <- keep_rows(d %>% filter(label == "manual"))
  H1 <- keep_rows(d %>% filter(label == "mmlb"))

  # (2) the mixed model's empirical probability limit, per (model, ntime),
  #     taken at the largest cluster count where Monte Carlo noise is smallest
  plim <- H1 %>% filter(method == "glmer") %>% group_by(model, ntime) %>%
    filter(N == max(N)) %>%
    summarise(plim = mean(estimate), n_plim = n(),
              mcse_plim = sd(estimate) / sqrt(n()), .groups = "drop") %>%
    mutate(dgp = lab)
  PL[[length(PL) + 1L]] <- plim

  # (3) intersection: replicates where every method returned a usable result.
  #     gee_ss contributes THREE rows (md, kc, fg), so count distinct methods.
  nmeth <- n_distinct(c(WALD, LRT))
  common <- H0 %>% distinct(model, N, ntime, sim_id, method) %>%
    count(model, N, ntime, sim_id, name = "nm") %>%
    filter(nm == nmeth) %>% select(model, N, ntime, sim_id)
  commonH1 <- H1 %>% distinct(model, N, ntime, sim_id, method) %>%
    count(model, N, ntime, sim_id, name = "nm") %>%
    filter(nm == nmeth) %>% select(model, N, ntime, sim_id)

  k <- k + 1L
  CW[[k]] <- bind_rows(
    H0 %>% group_by(method, model, N, ntime) %>%
      summarise(n = n(), k_rej = sum(p < 0.05), .groups = "drop") %>% mutate(set = "own"),
    H0 %>% inner_join(common, by = c("model", "N", "ntime", "sim_id")) %>%
      group_by(method, model, N, ntime) %>%
      summarise(n = n(), k_rej = sum(p < 0.05), .groups = "drop") %>% mutate(set = "intersection")
  ) %>% mutate(dgp = lab)

  # truth table built once per cell, then JOINED -- a rowwise lookup returning
  # NULL on any cell aborts the whole mutate, and it is slower besides.
  tt <- expand.grid(model = sort(unique(H1$model)), ntime = sort(unique(H1$ntime)))
  tt$t_cond <- as.numeric(mapply(cond_gen, lab, tt$model, tt$ntime))
  tt$t_marg <- as.numeric(mapply(marg_gen, lab, tt$model, tt$ntime))
  stopifnot(all(is.finite(tt$t_cond)), all(is.finite(tt$t_marg)))

  covtab <- function(df, setname) {
    df %>% filter(method %in% WALD) %>%
      left_join(plim %>% select(model, ntime, plim), by = c("model", "ntime")) %>%
      left_join(tt, by = c("model", "ntime")) %>%
      mutate(t_gen  = ifelse(method == "glmer", t_cond, t_marg),
             t_plim = ifelse(method == "glmer", plim,   t_marg)) %>%
      filter(is.finite(t_gen), is.finite(t_plim)) %>%
      mutate(c_gen  = abs(estimate - t_gen)  <= z * se,
             c_plim = abs(estimate - t_plim) <= z * se) %>%
      group_by(method, model, N, ntime) %>%
      summarise(n = n(),
                cov_gen  = 100 * mean(c_gen),
                cov_plim = 100 * mean(c_plim),
                mcse     = 100 * sqrt(mean(c_gen) * (1 - mean(c_gen)) / n()),
                .groups = "drop") %>% mutate(set = setname)
  }
  CV[[k]] <- bind_rows(
    covtab(H1, "own"),
    covtab(H1 %>% inner_join(commonH1, by = c("model", "N", "ntime", "sim_id")), "intersection")
  ) %>% mutate(dgp = lab)

  # empty-intersection cells, reported rather than averaged away
  allcells <- H0 %>% distinct(model, N, ntime)
  filled   <- common %>% distinct(model, N, ntime)
  SC[[k]]  <- allcells %>% anti_join(filled, by = c("model", "N", "ntime")) %>% mutate(dgp = lab)
  rm(d, H0, H1); gc(verbose = FALSE)
}
L <- bind_rows(CW); C <- bind_rows(CV); P <- bind_rows(PL); E <- bind_rows(SC)
write.csv(L, file.path(OUT, "H_level.csv"),         row.names = FALSE)
write.csv(C, file.path(OUT, "H_coverage.csv"),      row.names = FALSE)
write.csv(P, file.path(OUT, "H_glmer_plim.csv"),    row.names = FALSE)
write.csv(E, file.path(OUT, "H_empty_cells.csv"),   row.names = FALSE)

cat("===== (2) the mixed model's target vs where it actually converges =====\n")
print(as.data.frame(P %>% group_by(dgp) %>%
  summarise(plim_min = round(min(plim), 4), plim_max = round(max(plim), 4),
            .groups = "drop")), row.names = FALSE)

cat("\n===== COVERAGE, pooled, both targets and both denominators =====\n")
p2 <- C %>% group_by(dgp, method, set) %>%
  summarise(gen = sum(cov_gen * n) / sum(n), plim = sum(cov_plim * n) / sum(n),
            n = sum(n), .groups = "drop") %>%
  mutate(across(c(gen, plim), ~round(.x, 2)))
print(as.data.frame(p2 %>% filter(set == "own") %>% arrange(dgp, desc(gen))), row.names = FALSE)

cat("\n--- the mixed model only, both targets, by cluster count (T >= 4) ---\n")
g <- C %>% filter(method == "glmer", set == "own", ntime >= 4) %>%
  group_by(dgp, N) %>%
  summarise(gen = round(sum(cov_gen * n) / sum(n), 1),
            plim = round(sum(cov_plim * n) / sum(n), 1), .groups = "drop")
print(as.data.frame(g %>% pivot_wider(names_from = N, values_from = c(gen, plim))),
      row.names = FALSE)

cat("\n===== (3) cells with an EMPTY intersection =====\n")
cat(sprintf("%d of %d design cells\n", nrow(E), n_distinct(L$model, L$N, L$ntime) * 3))
print(as.data.frame(E %>% count(dgp, model, ntime, name = "cells")), row.names = FALSE)

cat("\n===== TYPE I ERROR, pooled, both denominators =====\n")
p1 <- L %>% group_by(dgp, method, set) %>%
  summarise(rate = 100 * sum(k_rej) / sum(n), n = sum(n), .groups = "drop") %>%
  mutate(rate = round(rate, 2)) %>%
  pivot_wider(names_from = set, values_from = c(rate, n))
print(as.data.frame(p1 %>% arrange(dgp, rate_own)), row.names = FALSE)

cat("\nwrote H_level.csv, H_coverage.csv, H_glmer_plim.csv, H_empty_cells.csv\n")
