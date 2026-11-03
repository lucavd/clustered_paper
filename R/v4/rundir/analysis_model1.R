# R/analysis_model1.R

analyze_model1_mmm <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
    mod <- mm(
      binY ~ time*trt,
      lv.formula = ~ 1,
      t.formula  = ~ 1,
      data = data,
      id   = id,
      q    = 20
    )
    s <- summary(mod)$mean.table
    data.frame(
      term     = rownames(s),
      estimate = s[,1],
      se       = s[,2],
      stat     = s[,3],
      p        = s[,4],
      method   = "mmm",
      N        = N,
      ntime    = ntime,
      label    = label
    )
 }

analyze_model1_gee <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  data$id <- factor(data$id)
    mod <- geeglm(
      binY ~ time*trt,
      id = id,
      data = data,
      family = binomial,
      corstr = "exchangeable",
      scale.fix = TRUE
    )
    s <- summary(mod)$coefficients
    data.frame(
      term     = rownames(s),
      estimate = s[,1],
      se       = s[,2],
      stat     = s[,3],
      p        = s[,4],
      method   = "gee",
      N        = N,
      ntime    = ntime,
      label    = label
    )
}

analyze_model1_glmer <- function(data, N, ntime, label) {
  data <- as.data.frame(data)
  data$id <- factor(data$id)
    mod <- glmer(
      binY ~ time*trt + (1 | id),
      data = data,
      family = binomial
    )
    s <- summary(mod)$coefficients
    D <- VarCorr(mod)$id[[1]]
    data.frame(
      term     = rownames(s),
      estimate = s[,1],
      se       = s[,2],
      stat     = s[,3],
      p        = s[,4],
      method   = "glmer",
      D        = D,
      N        = N,
      ntime    = ntime,
      label    = label
    )
}

# Shrinkage-augmented GLM: rescale the conditional (GLMM) coefficient onto the
# marginal scale by beta_m ~= beta_c / sqrt(1 + c^2 * D), c = 16*sqrt(3)/(15*pi).
#
# NB: this previously rescaled `estimate` ONLY, leaving `se`, `stat` and `p`
# copied verbatim from the glmer fit. The shrinkage method therefore had no
# inference of its own -- its p-values WERE the GLMM's, by construction, which is
# why its type I error was identical to the GLMM's to the decimal (11.0/11.5/...)
# across all scenarios. Any conclusion drawn from that was a tautology of the
# code, not a result. Fixed 2026-08-05: the SE is now propagated by the delta
# method, treating D as fixed, and stat/p are recomputed.
#
# IMPORTANT consequence, worth stating explicitly in the paper: because estimate
# and se are rescaled by the SAME factor, the Wald statistic estimate/se is
# scale-invariant, so the p-value is UNCHANGED by the shrinkage correction.
# The shrinkage GLM therefore has exactly the GLMM's type I error -- but now for
# a mathematical reason (invariance), not because of a coding bug. What DOES
# change is the confidence interval: it is now centred on the marginal estimand
# with a correspondingly rescaled width, so CI coverage against beta_m becomes
# meaningful and must be reported (referee comment 4 asked precisely for this).
#
# Caveat to report: this first-order delta method treats D-hat as fixed and
# ignores its sampling variability. With few clusters that is not negligible; a
# parametric bootstrap over D would be the fuller treatment.
analyze_model1_shrink <- function(glmer_results) {

  D <- glmer_results$D[1]
  c_const    <- 16 * sqrt(3) / (15 * pi)
  shrink     <- sqrt(1 + c_const^2 * D)
  inv_shrink <- 1 / shrink

  glmer_results |>
    mutate(
      estimate = estimate * inv_shrink,
      se       = se * inv_shrink,          # delta method, D treated as fixed
      stat     = estimate / se,            # unchanged by the common rescaling
      p        = 2 * pnorm(-abs(stat)),
      method = "shrink",
      shrink = shrink,
      inv_shrink = inv_shrink
    )
}


# summarize_results_model1 <- function(results) {
#   results |>
#     group_by(N, ntime, method, term, label) |>
#     summarise(
#       mean_est    = mean(estimate, na.rm = TRUE),
#       mean_se     = mean(se, na.rm = TRUE),
#       convergence = mean(!is.na(estimate)),
#       power       = mean(p < 0.05, na.rm = TRUE),
#       .groups     = "drop"
#     )
# }
