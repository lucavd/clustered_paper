# =============================================================================
# s7d_probe.R <scen> <dgp> -- S7d preparation (before pre-registration; organisation only, no measure enters the paper).
# (1) doubt (a) of S7c: does geesmv with corstr = "AR-M" return a genuine correction? Compares geesmv md/kc/fg with
#     fewclust 0.2.0 (laptop working tree, NOT a registered commit) and clubSandwich CR2/CR3 on the geeglm fit of the
#     same working correlation, and records the gee::gee refit that geesmv performs internally (beta, alpha).
# (2) serial per-dataset timing of every operation of the planned S7d run, data generation included, in each scenario.
# Points: N in {10, 25, 145} (M1c: 10, 25, 45), T in {2, 3, 10} (M1c: 2, 3, 6), i = 1..5, H0 ("manual").
# One process per (scen, dgp), serial inside: timings are wall-clock per operation.
# =============================================================================
a <- commandArgs(TRUE); SCEN <- a[1]; DGP <- a[2]
source(path.expand("~/MMM/R/s3/s3_common.R"))
source(path.expand("~/MMM/R/s7c/s7c_tests.R")); source(path.expand("~/MMM/R/s7d/s7d_flip.R"))
suppressMessages(library(data.table)); data.table::setDTthreads(1L)
.libPaths(c(.libPaths(), path.expand("~/R/s7d_probe_lib"), path.expand("~/R/s7c_lib")))
suppressMessages({ library(clubSandwich); library(fewclust); library(gee) })
env <- load_env(DGP, RUNDIR_V4)
if (SCEN %in% c("M6", "M7")) {
  source(path.expand("~/MMM/R/s6b/s6b_ns_common.R"))
  sl <- NS_SCEN[dgp == DGP & k == as.integer(sub("M", "", SCEN))]$slope; GEN <- make_gen(env, DGP, sl)
}
model <- if (SCEN == "M1c") 1L else as.integer(sub("M", "", SCEN))
FORM <- switch(SCEN, M2 = binY ~ time * trt + age, M5 = binY ~ time * trt + dose, binY ~ time * trt)
NS <- if (SCEN == "M1c") c(10L, 25L, 45L) else c(10L, 25L, 145L); TS <- if (SCEN == "M1c") c(2L, 3L, 6L) else c(2L, 3L, 10L)
CORS <- c(exch = "exchangeable", ind = "independence", ar1 = "ar1"); SMV <- c(exch = "exchangeable", ind = "independence", ar1 = "AR-M")
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err")); list(v = v, t = proc.time()[["elapsed"]] - t0) }
bad <- function(v) is.null(v) || inherits(v, "err") || inherits(v, "hard_timeout") || inherits(v, "fit_err")
msg <- function(v) if (inherits(v, "hard_timeout")) "timeout" else if (inherits(v, c("err", "fit_err"))) substr(as.character(v), 1, 100) else NA_character_
rows <- list(); add <- function(...) rows[[length(rows) + 1L]] <<- data.table(...)
for (N in NS) for (TT in TS) for (i in 1:5) {
  key <- list(scen = SCEN, dgp = DGP, N = N, ntime = TT, i = i)
  set_rep_seed(model, i)
  g <- tm(if (SCEN %in% c("M6", "M7")) GEN(N, TT, i, "manual") else gen_data(env, SCEN, model, N, TT, i, "manual"))
  do.call(add, c(key, list(op = "generate", cor = NA, sec = g$t, ok = !bad(g$v), se = NA_real_, est = NA_real_, alpha = NA_real_, msg = msg(g$v))))
  if (bad(g$v)) next
  d <- as.data.frame(g$v); d <- d[order(d$id, d$time), ]; key$n_obs <- nrow(d)
  for (cn in names(CORS)) {
    f <- tm(geepack::geeglm(FORM, id = id, data = d, family = binomial, corstr = CORS[[cn]], scale.fix = TRUE))
    fit <- if (bad(f$v)) NULL else f$v; k <- if (!is.null(fit)) match("time:trt", names(coef(fit))) else NA
    al <- if (!is.null(fit) && cn != "ind") as.numeric(fit$geese$alpha[1]) else NA_real_
    do.call(add, c(key, list(op = "geeglm", cor = cn, sec = f$t, ok = !is.null(fit), se = if (!is.null(fit)) sqrt(vcov(fit)[k, k]) else NA_real_,
                             est = if (!is.null(fit)) coef(fit)[[k]] else NA_real_, alpha = al, msg = msg(f$v))))
    gf <- tm(env$.with_hard_timeout(tryCatch(gee::gee(FORM, id = id, data = d, family = binomial, corstr = SMV[[cn]], silent = TRUE),
                                             error = function(e) structure(conditionMessage(e), class = "fit_err")), secs = 5))
    gk <- if (!bad(gf$v)) match("time:trt", names(gf$v$coefficients)) else NA
    do.call(add, c(key, list(op = "gee_refit", cor = cn, sec = gf$t, ok = !bad(gf$v), se = if (!bad(gf$v)) sqrt(gf$v$robust.variance[gk, gk]) else NA_real_,
                             est = if (!bad(gf$v)) gf$v$coefficients[[gk]] else NA_real_, alpha = if (!bad(gf$v) && cn != "ind") gf$v$working.correlation[1, 2] else NA_real_, msg = msg(gf$v))))
    for (ty in c("md", "kc", "fg")) {
      s <- tm(env$.with_hard_timeout(tryCatch(do.call(getExportedValue("geesmv", paste0("GEE.var.", ty)), list(FORM, id = "id", family = binomial, data = d, corstr = SMV[[cn]])),
                                             error = function(e) structure(conditionMessage(e), class = "fit_err")), secs = 5))
      kk <- match("time:trt", colnames(model.matrix(FORM, d)))
      do.call(add, c(key, list(op = paste0("geesmv_", ty), cor = cn, sec = s$t, ok = !bad(s$v), se = if (!bad(s$v)) sqrt(s$v$cov.beta[kk]) else NA_real_, est = NA_real_, alpha = NA_real_, msg = msg(s$v))))
      if (!is.null(fit)) { v <- tm(fewclust::gee_vcov(fit, type = ty))
        do.call(add, c(key, list(op = paste0("fewclust_", ty), cor = cn, sec = v$t, ok = !bad(v$v), se = if (!bad(v$v)) sqrt(v$v[k, k]) else NA_real_, est = NA_real_, alpha = NA_real_, msg = msg(v$v)))) }
    }
    if (!is.null(fit)) for (ty in c("CR2", "CR3")) { v <- tm(clubSandwich::vcovCR(fit, cluster = d$id, type = ty))
      do.call(add, c(key, list(op = paste0("club_", ty), cor = cn, sec = v$t, ok = !bad(v$v), se = if (!bad(v$v)) sqrt(as.matrix(v$v)[k, k]) else NA_real_, est = NA_real_, alpha = NA_real_, msg = msg(v$v)))) }
  }
  set.seed(20261022L + i); fl <- tm(s7d_flip(d, SCEN))
  do.call(add, c(key, list(op = "flip", cor = "ind", sec = fl$t, ok = !bad(fl$v) && is.finite(fl$v$p[1]), se = NA_real_, est = if (!bad(fl$v)) fl$v$p[1] else NA_real_, alpha = NA_real_, msg = if (!bad(fl$v)) fl$v$why[1] else msg(fl$v))))
}
R <- rbindlist(rows, fill = TRUE)
out <- path.expand(sprintf("~/MMM/results_v4/S7d/probe/probe_%s_%s.csv", SCEN, gsub("[^A-Za-z0-9]", "", DGP)))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE); fwrite(R, out)
cat(SCEN, DGP, "rows", nrow(R), "total sec", round(sum(R$sec), 1), "\n")
