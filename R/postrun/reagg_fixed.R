# =============================================================================
# reagg_fixed.R -- re-aggregate the stored run with the three code-review fixes.
#
# A1  scenario-5 GEE rows come from results_m5gee/ (refitted with sex+steroids)
#     rather than from the production chunks, where the uncorrected GEE had been
#     routed to scenario 4's fitter.
# E9  the exclusion filter is applied on the quantities each method's inference
#     actually uses: the likelihood-ratio test never touches `se`.
# E11 two denominators are reported side by side -- each method's own retained
#     set, and the INTERSECTION set on which every method returned a usable
#     result -- so a method that fails on the hardest cells cannot look better
#     calibrated by being scored on an easier subset.
# D1/D2 the conditional truth is joined per (dgp, model, ntime) from the values
#     the generators actually calibrate to, instead of a hard-coded 0.232/0.058.
#
# No simulation is re-run: this reads the 28,800 stored chunks plus the 5,760
# scenario-5 GEE chunks.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr)})
RES  <- path.expand("~/MMM/results_full")
M5   <- path.expand("~/MMM/results_m5gee")
OUT  <- path.expand("~/MMM/tables")
DGPS <- c(GenBinaryY = "GenBinaryY", GLMMnative = "GLMM-native", GLMMAR1 = "GLMM-AR1")

WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink")
LRT  <- "glmer_lrt"

# ---- truth tables --------------------------------------------------------
# marginal target per scenario (unchanged: this is what the marginal methods target)
MARG <- c("1" = 0.200, "2" = 0.200, "3" = 0.200, "4" = 0.050, "5" = 0.200)

ar1 <- read.csv(file.path(OUT, "F_truth_ar1.csv"))
ar1_truth <- ar1 %>% group_by(model, ntime) %>%
  summarise(cond = mean(cond_used), .groups = "drop") %>% mutate(dgp = "GLMM-AR1")

gn <- read.csv(file.path(OUT, "F_truth_glmmnative.csv"))   # implied MARGINAL, drifts with T
gn_marg <- gn %>% transmute(ntime, marg_main = m_main, marg_m4 = m_m4)

cond_truth <- function(dgp, model, ntime) {
  if (dgp == "GLMM-AR1") {
    v <- ar1_truth$cond[ar1_truth$model == model & ar1_truth$ntime == ntime]
    if (length(v)) return(v[1])
  }
  if (model == 4) 0.058 else 0.232          # GenBinaryY / GLMM-native, as generated
}
marg_truth <- function(dgp, model, ntime) {
  if (dgp == "GLMM-native") {
    r <- gn_marg[gn_marg$ntime == ntime, ]
    if (nrow(r)) return(if (model == 4) r$marg_m4[1] else r$marg_main[1])
  }
  unname(MARG[as.character(model)])
}

keep_rows <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method == LRT, converged, is.finite(p))
)

lvl <- list(); cvg <- list(); k <- 0L
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs <- list.files(RES, pattern = paste0("^", pref, "_"), full.names = TRUE)
  key <- sub("_c[0-9]+[.]rds$", "", basename(fs))
  full <- names(which(table(key) == 10))
  d <- bind_rows(lapply(fs[key %in% full], readRDS))
  d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt")

  # --- A1: swap in the refitted scenario-5 uncorrected GEE -----------------
  m5f <- list.files(M5, pattern = paste0("^", pref, "_"), full.names = TRUE)
  if (length(m5f)) {
    m5 <- bind_rows(lapply(m5f, readRDS)) %>% filter(term == "time:trt")
    m5$model <- 5L
    n_before <- sum(d$model == 5L & d$method == "gee")
    d <- d %>% filter(!(model == 5L & method == "gee"))
    d <- bind_rows(d, m5 %>% select(any_of(names(d))))
    cat(sprintf("%-12s scenario-5 GEE rows replaced: %d -> %d\n", lab, n_before,
                sum(d$model == 5L & d$method == "gee")))
  }

  H0 <- keep_rows(d %>% filter(label == "manual"))
  H1 <- keep_rows(d %>% filter(label == "mmlb"))

  # --- E11: each method's own set, AND the intersection set ---------------
  common <- H0 %>% count(model, N, ntime, sim_id, name = "nm") %>%
    filter(nm == length(c(WALD, LRT))) %>% select(model, N, ntime, sim_id)
  k <- k + 1L
  lvl[[k]] <- bind_rows(
    H0 %>% group_by(method, N) %>%
      summarise(n = n(), rate = 100 * mean(p < 0.05),
                mcse = 100 * sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / n()),
                .groups = "drop") %>% mutate(set = "own"),
    H0 %>% inner_join(common, by = c("model", "N", "ntime", "sim_id")) %>%
      group_by(method, N) %>%
      summarise(n = n(), rate = 100 * mean(p < 0.05),
                mcse = 100 * sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / n()),
                .groups = "drop") %>% mutate(set = "intersection")
  ) %>% mutate(dgp = lab)

  # --- D1/D2: coverage against the truth each mechanism actually used ------
  cvg[[k]] <- H1 %>% filter(method %in% WALD) %>%
    rowwise() %>%
    mutate(truth = if (method == "glmer") cond_truth(lab, model, ntime)
                   else marg_truth(lab, model, ntime)) %>%
    ungroup() %>%
    mutate(covers = abs(estimate - truth) <= 1.959964 * se) %>%
    group_by(method, N) %>%
    summarise(n = n(), cvg = 100 * mean(covers),
              mcse = 100 * sqrt(mean(covers) * (1 - mean(covers)) / n()),
              .groups = "drop") %>% mutate(dgp = lab)
  rm(d, H0, H1); gc(verbose = FALSE)
}
L <- bind_rows(lvl); C <- bind_rows(cvg)
write.csv(L, file.path(OUT, "G_level_fixed.csv"), row.names = FALSE)
write.csv(C, file.path(OUT, "G_coverage_fixed.csv"), row.names = FALSE)

cat("\n===== TYPE I ERROR, pooled (nominal 5%) =====\n")
cat("'own' = each method's own retained set; 'intersection' = replicates where\n")
cat("every method returned a usable result.\n\n")
p1 <- L %>% group_by(dgp, method, set) %>%
  summarise(rate = sum(rate * n) / sum(n), n = sum(n), .groups = "drop") %>%
  mutate(rate = round(rate, 2)) %>%
  pivot_wider(names_from = set, values_from = c(rate, n))
print(as.data.frame(p1 %>% arrange(dgp, rate_own)), row.names = FALSE)

cat("\n===== COVERAGE with the corrected truth (nominal 95%) =====\n")
p2 <- C %>% group_by(dgp, method) %>%
  summarise(cvg = round(sum(cvg * n) / sum(n), 2), n = sum(n), .groups = "drop")
print(as.data.frame(p2 %>% arrange(dgp, desc(cvg))), row.names = FALSE)

cat("\n===== scenario-5 GEE: before vs after the dispatch fix =====\n")
cat("(the level below now comes from a GEE fitted with sex and steroids)\n")
cat("\nwrote G_level_fixed.csv, G_coverage_fixed.csv\n")
