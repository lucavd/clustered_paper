# =============================================================================
# performance_metrics.R
# ADEMP performance measures with Monte Carlo standard errors.
#
# Follows Morris, White & Crowther (2019, Statistics in Medicine 38:2074-2102),
# "Using simulation studies to evaluate statistical methods". Every measure is
# reported with its Monte Carlo SE, which the submitted version omitted -- and
# which matters for how the results may be read: with 1000 replicates the MCSE
# on a proportion near 0.05 is ~0.7 percentage points, so the previously
# highlighted contrast "MMM 5.4% vs GEE 5.9%" is well inside noise and cannot
# be interpreted as a difference, whereas "GLMM 11.0% vs 5%" is far outside it.
#
# The review also asked for bias, empirical SD and RMSE alongside MAPE, noting
# that MAPE destabilises when the true value is small -- visible in the >200%
# values for M4, where beta = 0.05. MAPE is retained for continuity with the
# submitted tables but is no longer the primary accuracy measure.
#
# ESTIMAND HANDLING: methods target different estimands. Each row must carry
# the true value appropriate to ITS OWN estimand (`beta_true`), and where a
# common marginal comparison is wanted, `beta_true_marginal` as well. The
# review's first comment was precisely that scoring the GLMM against the
# marginal truth conflates estimand mismatch with estimation error.
# =============================================================================


#' ADEMP performance measures for one group of replicates
#'
#' @param est   vector of point estimates across replicates
#' @param se    vector of model-based standard errors
#' @param p     vector of p-values
#' @param truth true value of the estimand THIS method targets
#' @param alpha nominal significance level
#' @return one-row data frame of measures, each with its Monte Carlo SE
#' @param max_abs exclusion threshold for separated fits (see below)
#' @param uses_se whether this method's inference depends on the standard error.
#'   FALSE for a likelihood-ratio test, whose p-value comes from a deviance
#'   difference and never touches `se`. Filtering such rows on `se` discards
#'   results that are perfectly valid, and the discarded ones are not a random
#'   subset: separated replicates are disproportionately those that reject, so
#'   the exclusion understates the method's rejection rate. Measured on the
#'   stored run, the uniform filter dropped 113,791 of 1,435,382 finite
#'   likelihood-ratio rows (7.93%) -- the largest exclusion of any method.
ademp_measures <- function(est, se, p, truth, alpha = 0.05, max_abs = 10,
                           uses_se = TRUE) {
  finite <- is.finite(est) &
            (if (uses_se) is.finite(se) & se > 0 else TRUE)
  n_finite <- sum(finite)

  # SEPARATION EXCLUSION -- explicit, quantified, and reported.
  #
  # On sparse cells the outcome is completely separated by the covariates and the
  # MLE diverges: measured on the first 60 completed cells, |beta-hat| exceeded
  # 1e6 in 5.1% of uncorrected GEE fits, with a maximum of 5.9e+41 and SEs up to
  # 4e+30. A single such replicate destroys every moment-based measure: the mean
  # GEE estimate over all fits was -1.0e+37 against a true 0.2, while the mean
  # over non-separated fits was 0.097 and the median 0.096.
  #
  # The submitted analysis applied the same |estimate| < 10 filter (see
  # ci_coverage_analysis.R) but never reported it. Here the excluded fraction is
  # returned as `pct_separated` for every cell so it appears in the tables: it is
  # not a nuisance to hide but a property of the small-cluster regime the paper
  # is about (16.9% at N=10 falling to 1.8% at N=45).
  #
  # Exclusion is applied on the quantities each method's inference actually
  # uses -- uniform in EFFECT, not merely in form. An earlier version applied
  # the `se` conditions to every method, which was uniform in form and biased
  # in effect (see `uses_se`).
  keep <- finite & abs(est) < max_abs &
          (if (uses_se) se < max_abs else TRUE)
  n_excl <- n_finite - sum(keep)

  est <- est[keep]; se <- se[keep]; p <- p[keep]
  n <- length(est)
  if (n < 2L) {
    return(data.frame(n_rep = n, n_p = sum(!is.na(p)), n_finite = n_finite, n_separated = n_excl,
                      pct_separated = if (n_finite > 0) 100 * n_excl / n_finite else NA_real_,
                      bias = NA_real_, bias_mcse = NA_real_,
                      emp_se = NA_real_, emp_se_mcse = NA_real_,
                      model_se = NA_real_, se_ratio = NA_real_,
                      rmse = NA_real_, rmse_mcse = NA_real_,
                      mape = NA_real_, mape_mcse = NA_real_,
                      coverage = NA_real_, coverage_mcse = NA_real_,
                      rejection = NA_real_, rejection_mcse = NA_real_))
  }

  mean_est <- mean(est)
  bias     <- mean_est - truth
  # MCSE of bias is the empirical SD over sqrt(n) (Morris et al., Table 6)
  emp_se   <- sd(est)
  bias_mcse <- emp_se / sqrt(n)
  # MCSE of the empirical SE
  emp_se_mcse <- emp_se / sqrt(2 * (n - 1))

  model_se <- sqrt(mean(se^2))          # root mean model-based variance
  se_ratio <- model_se / emp_se         # <1 means SEs are too small

  rmse      <- sqrt(mean((est - truth)^2))
  rmse_mcse <- sd((est - truth)^2) / (2 * rmse * sqrt(n))

  # MAPE retained for continuity; unstable for |truth| near zero
  ape       <- if (abs(truth) > 1e-8) 100 * abs((est - truth) / truth) else rep(NA_real_, n)
  mape      <- mean(ape)
  mape_mcse <- sd(ape) / sqrt(n)

  z        <- qnorm(1 - alpha / 2)
  covered  <- (truth >= est - z * se) & (truth <= est + z * se)
  coverage <- mean(covered)
  coverage_mcse <- sqrt(coverage * (1 - coverage) / n)

  # The rejection rate is a mean over the replicates that HAVE a p-value, which
  # need not be all `n` retained rows; its MCSE must use that same denominator.
  n_p <- sum(!is.na(p))
  rejection <- mean(p < alpha, na.rm = TRUE)
  rejection_mcse <- if (n_p > 0) sqrt(rejection * (1 - rejection) / n_p) else NA_real_

  data.frame(
    n_rep = n,
    n_p = n_p,
    n_finite = n_finite,
    n_separated = n_excl,
    pct_separated = if (n_finite > 0) 100 * n_excl / n_finite else NA_real_,
    bias = bias, bias_mcse = bias_mcse,
    emp_se = emp_se, emp_se_mcse = emp_se_mcse,
    model_se = model_se, se_ratio = se_ratio,
    rmse = rmse, rmse_mcse = rmse_mcse,
    mape = mape, mape_mcse = mape_mcse,
    coverage = coverage, coverage_mcse = coverage_mcse,
    rejection = rejection, rejection_mcse = rejection_mcse
  )
}


#' Estimand lookup: the true value each method targets, per DGP and scenario
#'
#' This is the table the review asked to see stated explicitly ("a brief table
#' stating the target estimand and numerical true value for each method under
#' each data-generating process"). Marginal methods (MMM, GEE, shrinkage GLM)
#' target beta_m; the GLMM targets beta_c.
#'
#' beta_c values are the conditional parameters actually used by each DGP's
#' generators; for the AR(1) DGP they are solved per replicate by
#' calibrate_beta_c_ar1(), so they are supplied at analysis time rather than
#' hard-coded here.
estimand_table <- function() {
  marg <- c("1" = 0.2, "2" = 0.2, "3" = 0.2, "4" = 0.05, "5" = 0.2)
  cond <- c("1" = 0.232, "2" = 0.232, "3" = 0.232, "4" = 0.058, "5" = 0.232)
  rows <- list()
  for (dgp in c("GenBinaryY", "GLMM-native", "GLMM-AR1")) {
    for (m in names(marg)) {
      for (meth in c("mmm", "gee", "shrink", "gee_md", "gee_kc", "gee_fg")) {
        rows[[length(rows) + 1]] <- data.frame(
          dgp = dgp, model = as.integer(m), method = meth,
          estimand = "marginal", beta_true = unname(marg[m]),
          beta_true_marginal = unname(marg[m]))
      }
      for (meth in c("glmer", "glmer_lrt")) {
        rows[[length(rows) + 1]] <- data.frame(
          dgp = dgp, model = as.integer(m), method = meth,
          estimand = "conditional", beta_true = unname(cond[m]),
          beta_true_marginal = unname(marg[m]))
      }
    }
  }
  do.call(rbind, rows)
}


#' Summarise a raw replicate-level results table into ADEMP measures
#'
#' @param results replicate-level data frame (term, estimate, se, p, method,
#'                N, ntime, label, model, dgp)
#' @param term_of_interest coefficient to summarise
#' @param scoring "own" scores each method against its own estimand;
#'                "marginal" scores every method against the marginal truth
#'                (the two together answer the review's first comment)
summarise_ademp <- function(results, term_of_interest = "time:trt",
                            scoring = c("own", "marginal"), alpha = 0.05) {
  scoring <- match.arg(scoring)
  et <- estimand_table()
  truth_col <- if (scoring == "own") "beta_true" else "beta_true_marginal"

  d <- results[results$term %in% term_of_interest & results$label == "mmlb", , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- merge(d, et[, c("dgp", "model", "method", "estimand", truth_col)],
             by = c("dgp", "model", "method"), all.x = TRUE)
  names(d)[names(d) == truth_col] <- "truth"

  keys <- unique(d[, c("dgp", "model", "method", "N", "ntime")])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    k <- keys[i, ]
    s <- d[d$dgp == k$dgp & d$model == k$model & d$method == k$method &
             d$N == k$N & d$ntime == k$ntime, , drop = FALSE]
    if (!nrow(s) || all(is.na(s$truth))) next
    out[[i]] <- cbind(k, scoring = scoring,
                      estimand = s$estimand[1], truth = s$truth[1],
                      ademp_measures(s$estimate, s$se, s$p, s$truth[1], alpha))
  }
  do.call(rbind, Filter(Negate(is.null), out))
}


#' Convergence and exclusion accounting, per method and cell
#'
#' The review's minor comments asked for these; the submitted pipeline discarded
#' failures silently (`error = "null"` in _targets.R) so they could not be
#' reported at all.
summarise_convergence <- function(results, term_of_interest = "time:trt") {
  d <- results[results$term %in% term_of_interest | is.na(results$term), , drop = FALSE]
  keys <- unique(d[, c("dgp", "model", "method", "N", "ntime")])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    k <- keys[i, ]
    s <- d[d$dgp == k$dgp & d$model == k$model & d$method == k$method &
             d$N == k$N & d$ntime == k$ntime, , drop = FALSE]
    out[[i]] <- cbind(k, data.frame(
      n_attempted = nrow(s),
      pct_converged = 100 * mean(s$converged, na.rm = TRUE),
      pct_warned = 100 * mean(!is.na(s$status) & s$converged, na.rm = TRUE),
      pct_na_estimate = 100 * mean(is.na(s$estimate)),
      pct_extreme_gt10 = 100 * mean(abs(s$estimate) > 10, na.rm = TRUE)
    ))
  }
  do.call(rbind, out)
}
