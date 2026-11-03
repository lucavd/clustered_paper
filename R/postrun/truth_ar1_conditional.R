# =============================================================================
# truth_ar1c.R -- D1: recover the conditional truth the AR(1) generator actually
# used, by calling calibrate_beta_c_ar1() with exactly the arguments
# .ar1_generate() passes it, on designs regenerated from the production seeds.
#
# performance_metrics.R hard-codes cond = 0.232 (model 4: 0.058) for all three
# mechanisms, while stating in its own header that the AR(1) values are supplied
# at analysis time. Nothing supplies them.
#
# No model is fitted. Cheap.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr)})
MSM <- path.expand("~/MMM/MMM_sim_manual")
sys.source(file.path(MSM, "gen_ar1_data.R"), envir = globalenv())
sys.source(file.path(MSM, "R_ar1", "generation_models_ar1.R"), envir = globalenv())

# the per-scenario (formula, re.formula, D, design-extra) exactly as the
# generate_data_model* wrappers define them
# Read back from generate_data_model*() in generation_models_ar1.R -- note the
# mean model for scenarios 3-5 is ~ time * trt with NO covariate; the covariate
# enters only through the random-effects structure.
SPEC <- list(
  "1" = list(f = ~ time * trt,       re = ~ 1,             D = 1.0, ex = "none"),
  "2" = list(f = ~ time * trt + age, re = ~ 1,             D = 1.0, ex = "age"),
  "3" = list(f = ~ time * trt,       re = ~ 1 + age,
             D = matrix(c(1.0, 0.01, 0.01, 0.001), 2, byrow = TRUE), ex = "age"),
  "4" = list(f = ~ time * trt,       re = ~ 1 + age,
             D = matrix(c(1.0, 0.01, 0.01, 0.001), 2, byrow = TRUE), ex = "age"),
  "5" = list(f = ~ time * trt,       re = ~ 1 + age + sex,
             D = matrix(c(1.0, 0.01, 0.01, 0.01, 0.001, 0.0, 0.01, 0.0, 0.01),
                        3, byrow = TRUE), ex = "full")
)

T_LEVELS <- c(2, 3, 4, 6, 8, 10)
N_LEVELS <- c(25, 75, 145)

rows <- list(); k <- 0L
for (m in as.character(1:5)) {
  sp <- SPEC[[m]]
  for (T in T_LEVELS) for (N in N_LEVELS) {
    vals <- numeric(0)
    for (rep in 1:6) {
      set.seed(20260805L + 1000L * as.integer(m) + rep)
      d <- tryCatch(.ar1_design(N, T, sp$ex), error = function(e) NULL)
      if (is.null(d)) next
      tgt <- tryCatch(.ar1_target(as.integer(m), "mmlb"), error = function(e) NULL)
      if (is.null(tgt)) next
      bc <- tryCatch(calibrate_beta_c_ar1(mean.formula = sp$f, re.formula = sp$re,
                                          target_beta_m = tgt, D = sp$D,
                                          sigma2_ar = AR1_SIGMA2, data = d),
                     error = function(e) NULL)
      if (!is.null(bc)) {
        b <- as.numeric(bc)
        # interaction is the term named time:trt in this design matrix
        nm <- colnames(model.matrix(sp$f, model.frame(sp$f, d)))
        j  <- match("time:trt", nm)
        if (!is.na(j)) vals <- c(vals, b[j])
      }
    }
    if (length(vals)) {
      k <- k + 1L
      rows[[k]] <- data.frame(model = as.integer(m), ntime = T, N = N,
                              cond_used = mean(vals), sd_ = sd(vals), n_ok = length(vals))
    }
  }
}
A <- bind_rows(rows)
if (!nrow(A)) { cat("recovery failed -- inspect the calibration arguments\n"); quit(status = 0) }

cat("=== conditional time:trt the AR(1) generator actually calibrates to ===\n")
print(as.data.frame(A %>% group_by(model, ntime) %>%
  summarise(cond = round(mean(cond_used), 4), .groups = "drop") %>%
  pivot_wider(names_from = ntime, values_from = cond, names_prefix = "T=")),
  row.names = FALSE)

A$hardcoded <- ifelse(A$model == 4, 0.058, 0.232)
A$pct_error <- 100 * (A$hardcoded - A$cond_used) / A$cond_used
cat("\n=== the metrics module hard-codes 0.232 (model 4: 0.058). Error, percent ===\n")
print(as.data.frame(A %>% group_by(model) %>%
  summarise(min = round(min(pct_error), 1), max = round(max(pct_error), 1),
            mean = round(mean(pct_error), 1), .groups = "drop")), row.names = FALSE)
write.csv(A, path.expand("~/MMM/tables/F_truth_ar1.csv"), row.names = FALSE)
cat("\nwrote F_truth_ar1.csv\n")
