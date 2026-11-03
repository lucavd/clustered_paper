# =============================================================================
# s7d_common.R -- S7d: run with alternative working correlations, independent corrections and the sign-flip test
# on all 1,896 design points (main grid M1-M5 x 3 mechanisms, M1c x 3 mechanisms, M6-M7 x 2 mechanisms), H0 and H1,
# 1,000 replicates per arm, chunks of 100. Datasets regenerated exactly as in the study (R/s3/s3_common.R::one_rep,
# R/s6b/s6b_ns_common.R::one_rep_ns): degenerate-design flag first, then seeded generation (20260805 + 1000 k + i,
# L'Ecuyer-CMRG). Per replicate and working correlation (exchangeable, independence, ar1):
#   geepack::geeglm (scale.fix = TRUE, the study's mean model) -> estimate, uncorrected sandwich SE, alpha
#   fewclust::gee_vcov md / kc / fg (fewclust 0.2.0, commit in ~/R/s7d_lib/fewclust_commit.txt, D-49)
#   clubSandwich CR2 (exchangeable and independence only; S7d technical decision on AR(1))
#   geesmv md / kc / fg under the study's 5 s hard timeout, replicates 1-50 only, exchangeable and independence (D-48)
# and once per replicate the sign-flip score test, primary and strict versions (D-46, D-47; R/s7d/s7d_flip.R),
# after set.seed(20261022 + i) as in S7c.
# =============================================================================
source(path.expand("~/MMM/R/s6b/s6b_ns_common.R"))   # sources R/s3/s3_common.R; make_gen, NS_SCEN
source(path.expand("~/MMM/R/s7c/s7c_tests.R")); source(path.expand("~/MMM/R/s7d/s7d_flip.R"))
suppressMessages(library(data.table)); data.table::setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
S7D_LIB <- path.expand("~/R/s7d_lib"); .libPaths(c(.libPaths(), S7D_LIB))
suppressMessages({ library(clubSandwich); library(fewclust) })
OUT_S7D <- file.path(REPO, "results_v4/S7d"); dir.create(file.path(OUT_S7D, "rows"), recursive = TRUE, showWarnings = FALSE)
FLIP_SEED <- function(i) 20261022L + as.integer(i)
SMV_MAX_I <- 50L
CORS <- c(exch = "exchangeable", ind = "independence", ar1 = "ar1"); SMV <- c(exch = "exchangeable", ind = "independence")
FORM_OF <- function(scen) switch(scen, M2 = binY ~ time * trt + age, M5 = binY ~ time * trt + dose, binY ~ time * trt)

s7d_tasks <- function() {
  a <- CJ(dgp = DGPS, scen = sprintf("M%d", 1:5), N = as.integer(N_LEVELS), ntime = as.integer(T_LEVELS))
  b <- CJ(dgp = DGPS, scen = "M1c", N = as.integer(N_M1C), ntime = as.integer(T_M1C))
  cc <- CJ(s = seq_len(nrow(NS_SCEN)), N = as.integer(N_LEVELS), ntime = as.integer(T_LEVELS))
  cc <- cbind(NS_SCEN[cc$s, .(dgp, scen = paste0("M", k))], cc[, .(N, ntime)])
  x <- rbind(a, b, cc)
  stopifnot(nrow(x) == 1896L)
  x <- x[, .(label = c("mmlb", "manual")), by = .(dgp, scen, N, ntime)][, .(chunk = 1:10), by = .(dgp, scen, N, ntime, label)]
  x[, model := fifelse(scen == "M1c", 1L, as.integer(sub("M", "", scen)))]
  x[, file := file.path(OUT_S7D, "rows", sprintf("%s_%s_N%d_T%d_%s_c%02d.rds", dgp_tag(dgp), tolower(scen), N, ntime, label, chunk))]
  x
}
s7d_gen <- function(env, t) {
  if (t$scen == "M1c") function(N, ntime, i, label) env$generate_data_model1c(N, ntime, i, label)
  else if (t$scen %in% c("M6", "M7")) make_gen(env, t$dgp, NS_SCEN[dgp == t$dgp & k == t$model]$slope)
  else function(N, ntime, i, label) env$make_all_data(N, ntime, i, t$model, label)
}
.bad <- function(v) is.null(v) || inherits(v, c("error", "hard_timeout", "fit_err"))
.msg <- function(v) if (inherits(v, "hard_timeout")) "geesmv hard timeout" else if (inherits(v, "error")) substr(conditionMessage(v), 1, 160) else if (inherits(v, "fit_err")) substr(as.character(v), 1, 160) else NA_character_
se_of <- function(V, k) { s <- tryCatch(sqrt(as.matrix(V)[k, k]), error = function(e) NA_real_, warning = function(w) NA_real_); if (is.finite(s)) s else NA_real_ }

one_cor <- function(env, d, form, cn, i) {
  r <- data.table(cor = cn, gee_ok = FALSE, est = NA_real_, se_unc = NA_real_, alpha = NA_real_, gee_msg = NA_character_,
                  se_md = NA_real_, se_kc = NA_real_, se_fg = NA_real_, se_cr2 = NA_real_,
                  msg_md = NA_character_, msg_kc = NA_character_, msg_fg = NA_character_, msg_cr2 = NA_character_,
                  smv_md = NA_real_, smv_kc = NA_real_, smv_fg = NA_real_, smv_msg_md = NA_character_, smv_msg_kc = NA_character_, smv_msg_fg = NA_character_)
  fit <- tryCatch(geepack::geeglm(form, id = id, data = d, family = binomial, corstr = CORS[[cn]], scale.fix = TRUE), error = function(e) e)
  if (.bad(fit)) { r$gee_msg <- .msg(fit); fit <- NULL }
  if (!is.null(fit)) {
    k <- match("time:trt", names(coef(fit)))
    r$gee_ok <- TRUE; r$est <- unname(coef(fit)[k]); r$se_unc <- se_of(vcov(fit), k)
    if (cn != "ind") r$alpha <- as.numeric(fit$geese$alpha[1])
    for (ty in c("md", "kc", "fg")) {
      v <- tryCatch(fewclust::gee_vcov(fit, type = ty), error = function(e) e)
      if (.bad(v)) set(r, 1L, paste0("msg_", ty), .msg(v)) else set(r, 1L, paste0("se_", ty), se_of(v, k))
    }
    if (cn != "ar1") { v <- tryCatch(clubSandwich::vcovCR(fit, cluster = d$id, type = "CR2"), error = function(e) e)
      if (.bad(v)) r$msg_cr2 <- .msg(v) else r$se_cr2 <- se_of(v, k) }
  }
  if (i <= SMV_MAX_I && cn %in% names(SMV)) {
    kk <- match("time:trt", colnames(model.matrix(form, d)))
    for (ty in c("md", "kc", "fg")) {
      v <- env$.with_hard_timeout(tryCatch(do.call(getExportedValue("geesmv", paste0("GEE.var.", ty)),
             list(form, id = "id", family = binomial, data = d, corstr = SMV[[cn]])),
             error = function(e) structure(conditionMessage(e), class = "fit_err")), secs = 5)
      if (.bad(v)) set(r, 1L, paste0("smv_msg_", ty), .msg(v)) else set(r, 1L, paste0("smv_", ty), sqrt(v$cov.beta[kk]))
    }
  }
  r
}
one_rep_s7d <- function(env, gen, t, i) {
  deg <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) e)
  rec <- data.table(sim_id = i, generated = !inherits(d, "error"), design_degenerate = deg,
                    gen_error = if (inherits(d, "error")) substr(conditionMessage(d), 1, 180) else NA_character_,
                    n_obs = NA_integer_, sum_y = NA_real_, p_flip = NA_real_, p_flip_strict = NA_real_, flip_why = NA_character_,
                    flip_strict_why = NA_character_, flip_max_abs_coef = NA_real_, flip_n_ref = NA_integer_)
  if (inherits(d, "error")) return(list(rec = rec, rows = NULL))
  d <- d[order(d$id, d$time), ]; rec$n_obs <- nrow(d); rec$sum_y <- sum(d$binY)
  form <- FORM_OF(t$scen)
  rows <- rbindlist(lapply(names(CORS), function(cn) one_cor(env, d, form, cn, i)))
  rows[, sim_id := i]
  set.seed(FLIP_SEED(i))
  fl <- tryCatch(s7d_flip(d, t$scen), error = function(e) NULL)
  if (is.null(fl)) rec$flip_why <- "flip error" else {
    rec$p_flip <- fl$p[1]; rec$p_flip_strict <- fl$p[2]; rec$flip_why <- fl$why[1]; rec$flip_strict_why <- fl$why[2]
    rec$flip_max_abs_coef <- fl$max_abs_coef[1]; rec$flip_n_ref <- fl$n_ref[1] }
  list(rec = rec, rows = rows)
}
run_chunk_s7d <- function(t) {
  env <- load_env(t$dgp, RUNDIR_V4); gen <- s7d_gen(env, t)
  ids <- (t$chunk - 1L) * 100L + seq_len(100L)
  out <- lapply(ids, function(i) one_rep_s7d(env, gen, t, i))
  recs <- rbindlist(lapply(out, `[[`, "rec"), fill = TRUE); rows <- rbindlist(lapply(out, `[[`, "rows"), fill = TRUE)
  for (D in list(recs, rows)) if (nrow(D)) D[, `:=`(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
  tmp <- paste0(t$file, ".tmp"); saveRDS(list(recs = recs, rows = rows), tmp); file.rename(tmp, t$file)
  nrow(rows)
}
