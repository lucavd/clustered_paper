# NB: scenario 5 previously routed to analyze_model4_* (MMM) and analyze_model3_*
# (GLMM), so under M5 the MMM was fitted without sex in the latent structure and
# without steroids in the transition model, and the GLMM without sex/steroids
# altogether -- three different specifications compared as if they were the same
# scenario. analyze_model5_* existed but was never called. Fixed 2026-08-05.
do_analysis_mmm <- function(data, N, ntime, label, model) {
  model |>
    switch(
      "1" = analyze_model1_mmm(data, N, ntime, label),
      "2" = analyze_model2_mmm(data, N, ntime, label),
      "3" = analyze_model3_mmm(data, N, ntime, label),
      "4" = analyze_model4_mmm(data, N, ntime, label),
      "5" = analyze_model5_mmm(data, N, ntime, label)
    )
}

do_analysis_gee <- function(data, N, ntime, label, model) {
  model |>
    switch(
      "1" = analyze_model1_gee(data, N, ntime, label),
      "2" = analyze_model2_gee(data, N, ntime, label),
      "3" = analyze_model3_gee(data, N, ntime, label),
      "4" = analyze_model4_gee(data, N, ntime, label),
      "5" = analyze_model5_gee(data, N, ntime, label)
    )
}

do_analysis_glmer <- function(data, N, ntime, label, model) {
  model |>
    switch(
      "1" = analyze_model1_glmer(data, N, ntime, label),
      "2" = analyze_model2_glmer(data, N, ntime, label),
      "3" = analyze_model3_glmer(data, N, ntime, label),
      "4" = analyze_model4_glmer(data, N, ntime, label),
      "5" = analyze_model5_glmer(data, N, ntime, label)
    )
}

# Fits are wrapped so that a non-converging model yields a row carrying its
# failure reason instead of vanishing. Previously `error = "null"` in _targets.R
# made failed fits disappear silently, so per-method convergence rates could not
# be reported at all -- referee minor comment asked for exactly these numbers.
# Added 2026-08-05.
do_analysis <- function(
  .data, n, ntime, label, model, technique, sim_id
) {
  fit <- function(expr) {
    warns <- character(0)
    out <- withCallingHandlers(
      tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "fit_error")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    list(out = out, warns = warns)
  }

  # Mean-model formula per scenario, needed by the small-sample competitors
  # (they refit rather than post-process, so they need the same specification
  # the standard GEE/GLMM analyses use for this scenario).
  mean_rhs <- switch(
    as.character(model),
    "1" = "time*trt",
    "2" = "time*trt + age",
    "3" = "time*trt",
    "4" = "time*trt",
    "5" = "time*trt + dose"
  )
  f_gee   <- as.formula(paste("binY ~", mean_rhs))
  f_glmer <- as.formula(paste("binY ~", mean_rhs, "+ (1|id)"))

  r <- technique |>
    switch(
      "mmm"       = fit(do_analysis_mmm(.data, n, ntime, label, model)),
      "gee"       = fit(do_analysis_gee(.data, n, ntime, label, model)),
      "glmer"     = fit(do_analysis_glmer(.data, n, ntime, label, model)),
      # Small-sample corrected competitors (see R/analysis_smallsample.R).
      # gee_ss returns three rows at once (Mancl-DeRouen, Kauermann-Carroll,
      # Fay-Graubard) since the three corrections disagree materially.
      "gee_ss"    = fit(analyze_gee_smallsample(.data, f_gee, n, ntime, label)),
      "glmer_lrt" = fit(analyze_glmer_lrt(.data, f_glmer, n, ntime, label))
    )

  # A small-sample helper returning NULL means the fit could not be produced.
  if (is.null(r$out)) {
    return(tibble::tibble(
      term = NA_character_, estimate = NA_real_, se = NA_real_,
      stat = NA_real_, p = NA_real_, method = technique,
      N = n, ntime = ntime, label = label, model = model, sim_id = sim_id,
      converged = FALSE, status = "fit returned NULL"
    ))
  }

  if (inherits(r$out, "fit_error")) {
    return(tibble::tibble(
      term = NA_character_, estimate = NA_real_, se = NA_real_,
      stat = NA_real_, p = NA_real_, method = technique,
      N = n, ntime = ntime, label = label, model = model, sim_id = sim_id,
      converged = FALSE,
      status = substr(as.character(r$out), 1, 200)
    ))
  }

  # The small-sample helpers set `converged`/`status` themselves (a geesmv
  # timeout is a non-convergence for that correction). Do NOT clobber those --
  # blanket-assigning converged = TRUE here would have silently reported timed-out
  # fits as clean successes.
  warn_txt <- if (length(r$warns)) substr(paste(r$warns, collapse = " | "), 1, 200) else NA_character_
  aux <- r$out |> mutate(model = model)
  if (!"converged" %in% names(aux)) aux$converged <- TRUE
  if (!"status" %in% names(aux)) {
    # lme4 reports non-convergence as a warning on an object that still returns
    # coefficients; record it rather than treating the fit as clean.
    aux$status <- warn_txt
  } else {
    aux$status <- ifelse(is.na(aux$status), warn_txt, aux$status)
  }
  aux <- tibble::as_tibble(aux)

  res <- if (technique == "glmer") {
    dplyr::bind_rows(
      aux,
      analyze_model1_shrink(aux) |>
        mutate(model = model) |>
        tibble::as_tibble()
    )
  } else {
    aux
  }

  res |>
    mutate(sim_id = sim_id)
}
