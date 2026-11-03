# =============================================================================
# analysis_smallsample.R
# Small-sample corrected competitors, added in response to the BMC MRM review.
#
# Motivation: the submitted version compared MMM against a GEE using the
# uncorrected sandwich estimator and a GLMM using the Wald test. Both are known
# to behave poorly with few clusters, so concluding that "GEE is ill-suited to
# NICU studies" and that "GLMM inflates type I error" risked comparing MMM
# against deliberately weakened opponents. These functions add the standard
# remedies so the comparison is like-for-like.
#
# Verified available on the analysis machine (2026-08-05):
#   - geesmv::GEE.var.md  Mancl & DeRouen (2001) bias-corrected sandwich
#   - geesmv::GEE.var.kc  Kauermann & Carroll (2001)
#   - geesmv::GEE.var.fg  Fay & Graubard (2001)
#   - likelihood-ratio test for the GLMM (anova on nested fits)
# NOT usable here, checked and rejected:
#   - clubSandwich CR2: vcovCR has no method for glmerMod (S4); works on
#     glm/geeglm only
#   - pbkrtest::PBmodcomp: targets linear mixed models; errors on a binomial
#     glmer with "unused argument (REML = FALSE)"
# =============================================================================


#' Run an expression under a hard, process-level timeout
#'
#' WHY THIS IS NECESSARY. Measured on lesexp-server 2026-08-05: all three geesmv
#' correction routines enter an infinite loop on sparse data -- 8.0% of fits
#' overall (30/375), concentrated exactly in the new small-cluster cells:
#'
#'     N=10 T=2  ->  4/25 hang (datasets with 2-4 events)
#'     N=15 T=2  ->  4/25 hang (3-5 events)
#'     N=20 T=3  ->  2/25 hang (6-9 events)
#'     N=25 T=4  ->  0/25
#'     N=45 T=4  ->  0/25
#'
#' The loop is inside geesmv's COMPILED code, so `setTimeLimit()` never fires
#' (R time limits are only checked at R-level evaluation points) -- verified: a
#' run under setTimeLimit spun five workers at 99.9% CPU for 30 minutes and had
#' to be killed externally. Only forking the call and SIGKILLing the child
#' stops it. Without this guard the extended grid cannot run at all: one hung
#' fit permanently consumes a worker, and at an 8% rate on the sparse cells the
#' whole pipeline deadlocks.
#'
#' A timed-out fit is reported as a non-convergence for that method, which is
#' the honest accounting -- and the rate must be reported, since it is itself a
#' finding about the applicability of these corrections at small N.
.with_hard_timeout <- function(expr, secs = 30) {
  f <- function() expr
  p <- parallel::mcparallel(f())
  r <- parallel::mccollect(p, wait = FALSE, timeout = secs)
  if (is.null(r)) {
    tools::pskill(p$pid, tools::SIGKILL)
    parallel::mccollect(p, wait = FALSE)
    return(structure("timeout", class = "hard_timeout"))
  }
  r[[1]]
}


#' GEE with small-sample corrected sandwich variance
#'
#' Refits nothing: takes the point estimates from the standard GEE fit and
#' replaces the variance with a bias-corrected version, which is how these
#' corrections are meant to be used. Returns one row per correction so all
#' three can be compared, since they do NOT agree with each other (at N=60 we
#' measured md se=0.1534, kc se=0.1005, fg se=0.1552 on the same fit).
#'
#' Each correction is called under a hard timeout; see .with_hard_timeout().
#'
#' @param data  analysis data frame, must be sorted by id
#' @param formula model formula for the mean structure
#' @param term_of_interest coefficient to extract
#' @param timeout_secs per-correction wall-clock cap. 5s is ~10-70x the measured
#'   time of a SUCCESSFUL correction (0.07s at N=15,T=2; 0.21s at N=45,T=6;
#'   0.57s at N=145,T=10), so it cannot truncate a genuine fit -- while a hung
#'   one loops without bound. Profiled 2026-08-05: at 30s the sparse cells cost
#'   15.1 s/replicate (3 corrections x 30s worst case = 90s) and the full grid
#'   projected to ~26h; at 5s it projects to ~17h with identical results.
analyze_gee_smallsample <- function(data, formula, N, ntime, label,
                                    term_of_interest = "time:trt",
                                    timeout_secs = 5) {
  d <- as.data.frame(data)
  d <- d[order(d$id), ]

  base <- tryCatch(
    geepack::geeglm(formula, id = id, data = d, family = binomial,
                    corstr = "exchangeable", scale.fix = TRUE),
    error = function(e) NULL
  )
  if (is.null(base)) return(NULL)

  est_all <- coef(base)
  if (!term_of_interest %in% names(est_all)) return(NULL)
  k   <- match(term_of_interest, names(est_all))
  est <- unname(est_all[k])

  out <- list()
  for (nm in c("md", "kc", "fg")) {
    fn <- paste0("GEE.var.", nm)
    v <- .with_hard_timeout(
      tryCatch(
        do.call(getExportedValue("geesmv", fn),
                list(formula, id = "id", family = binomial, data = d,
                     corstr = "exchangeable")),
        error = function(e) structure(conditionMessage(e), class = "fit_err")
      ),
      secs = timeout_secs
    )

    # Record failures rather than dropping them: the timeout rate at small N is
    # itself a result about whether these corrections are usable there.
    if (inherits(v, "hard_timeout") || inherits(v, "fit_err") || is.null(v)) {
      out[[nm]] <- data.frame(
        term = term_of_interest, estimate = est, se = NA_real_,
        stat = NA_real_, p = NA_real_, method = paste0("gee_", nm),
        N = N, ntime = ntime, label = label, converged = FALSE,
        status = if (inherits(v, "hard_timeout")) "geesmv hard timeout"
                 else if (inherits(v, "fit_err")) substr(as.character(v), 1, 200)
                 else "geesmv returned NULL",
        stringsAsFactors = FALSE
      )
      next
    }

    se <- sqrt(v$cov.beta[k])
    if (!is.finite(se) || se <= 0) {
      out[[nm]] <- data.frame(
        term = term_of_interest, estimate = est, se = NA_real_,
        stat = NA_real_, p = NA_real_, method = paste0("gee_", nm),
        N = N, ntime = ntime, label = label, converged = FALSE,
        status = "non-finite or non-positive corrected variance",
        stringsAsFactors = FALSE
      )
      next
    }

    out[[nm]] <- data.frame(
      term = term_of_interest, estimate = est, se = se,
      stat = est / se, p = 2 * pnorm(-abs(est / se)),
      method = paste0("gee_", nm),
      N = N, ntime = ntime, label = label,
      converged = TRUE, status = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}


#' GLMM inference by likelihood-ratio test instead of Wald
#'
#' The Wald test is known to be anti-conservative with few clusters, so the
#' 11% type I error reported for GLMM in the submitted version may partly be a
#' property of the TEST rather than of the framework. This refits the model
#' without the tested term and compares by LRT, keeping the Wald point estimate
#' and standard error (the LRT supplies the p-value only).
analyze_glmer_lrt <- function(data, formula, N, ntime, label,
                              term_of_interest = "time:trt") {
  d <- as.data.frame(data)
  d$id <- factor(d$id)

  full <- tryCatch(suppressMessages(lme4::glmer(formula, data = d, family = binomial,
                     control = lme4::glmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 1e5)))),
                   error = function(e) NULL)
  if (is.null(full)) return(NULL)

  fe <- lme4::fixef(full)
  if (!term_of_interest %in% names(fe)) return(NULL)
  s  <- summary(full)$coefficients
  est <- s[term_of_interest, 1]; se <- s[term_of_interest, 2]

  # Build the reduced model by DROPPING THE TERM FROM THE FORMULA TEXT.
  #
  # Two earlier attempts were wrong, both worth recording:
  #
  # 1. update(full, . ~ . - time:trt) re-evaluates the stored call in the
  #    formula's environment. When the formula is built by as.formula(paste(...))
  #    in a caller (as do_analysis does), that environment is the caller's frame
  #    where the local `d` does not exist, so update() failed and the method
  #    returned NULL for EVERY replicate.
  #
  # 2. reformulate(term.labels) drops the parentheses around random-effect terms:
  #    terms() reports "(1|id)" as the bare label "1 | id", so reformulate
  #    produced `binY ~ time + trt + 1 | id` -- which lme4 reads as a random
  #    effect on (time + trt + 1), giving
  #      "number of observations (=64) < number of random effects (=120)".
  #    This failed on 99.9% of T=2 and 99.4% of T=3 cells (and ~4% at T=4),
  #    i.e. exactly the short-series regime the extended grid was added to probe,
  #    while looking fine at T>=6. A method that silently vanishes from the cells
  #    it was introduced to inform would have been the worst kind of gap.
  #
  # Correct approach: take terms(), drop the tested one, and rebuild -- putting
  # the parentheses back around any random-effect term. (Plain text substitution
  # does not work either: deparse gives "time * trt", not "time:trt", so a
  # pattern on the label matches nothing.)
  tl <- attr(terms(formula), "term.labels")
  if (!term_of_interest %in% tl) return(NULL)
  keep  <- tl[tl != term_of_interest]
  is_re <- grepl("\\|", keep)
  rhs   <- paste(c(keep[!is_re], paste0("(", keep[is_re], ")")), collapse = " + ")
  f_red <- as.formula(paste(all.vars(formula)[1], "~", rhs))
  environment(f_red) <- environment()

  red <- tryCatch(suppressMessages(lme4::glmer(f_red, data = d, family = binomial,
                    control = lme4::glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 1e5)))),
                  error = function(e) NULL)
  if (is.null(red)) return(NULL)

  a <- tryCatch(anova(red, full), error = function(e) NULL)
  if (is.null(a)) return(NULL)
  chisq <- a$Chisq[2]; p_lrt <- a$`Pr(>Chisq)`[2]
  if (!is.finite(p_lrt)) return(NULL)

  data.frame(
    term = term_of_interest, estimate = est, se = se,
    stat = chisq, p = p_lrt, method = "glmer_lrt",
    N = N, ntime = ntime, label = label,
    converged = TRUE, status = NA_character_,
    stringsAsFactors = FALSE
  )
}
