# =============================================================================
# s6b_ns_common.R -- S6b part B (D-34): scenarios with a large slope of the latent SD in time.
# Copies of M3 (mean binY ~ time * trt, latent SD with a slope in time, no transition covariate) in which ONLY the
# slope changes:
#   GenBinaryY  sigma = (1, 0.01) -> (1, 0.1) [M6] and (1, 0.25) [M7]          (latent SD 1 + s * t)
#   GLMM-native D = [1, 0.01; 0.01, 0.001] -> slope variance 0.01 [M6] and 0.04 [M7]
# The generator is the v4 generate_data_model3 of each family, with the slope substituted in its source text
# (checked: exactly the expected number of substitutions). Seeds 20260805 + 1000 * k + i, k = 6, 7 (no overlap with
# k = 1..5), L'Ecuyer-CMRG (D-07). One replicate as R/s3/s3_common.R::one_rep (degenerate-design flag first, then
# seeded generation, then the five analyses of the study with the M3 formulas), plus the MMM fits of the D-32 rule
# (fit A as in the study, fit B with lv.formula = ~1; threshold 1e-3).
# =============================================================================
source(path.expand("~/MMM/R/s3/s3_common.R"))   # libraries, load_env, set_rep_seed, design_degenerate, trace_row, N_LEVELS, T_LEVELS
suppressMessages(library(data.table))
OUT_NS <- file.path(REPO, "results_v4/S6b/ns")
NS_SCEN <- data.table(dgp = c("GenBinaryY", "GenBinaryY", "GLMM-native", "GLMM-native"), k = c(6L, 7L, 6L, 7L),
                      slope = c(0.1, 0.25, 0.01, 0.04),
                      what = c("sigma = (1, 0.1)", "sigma = (1, 0.25)", "D[2,2] = 0.01", "D[2,2] = 0.04"))
# source text substitution: returns the modified generate_data_model3 of env, evaluated in env
make_gen <- function(env, dgp, slope) {
  src <- deparse(env$generate_data_model3, width.cutoff = 500L)
  if (dgp == "GenBinaryY") { pat <- "c(1, 0.01)"; rep <- sprintf("c(1, %s)", format(slope)); n_exp <- 2L }
  else { pat <- "0.01, 0.001)"; rep <- sprintf("0.01, %s)", format(slope)); n_exp <- 1L }
  hits <- lengths(regmatches(src, gregexpr(pat, src, fixed = TRUE)))
  stopifnot(sum(hits) == n_exp)
  src <- gsub(pat, rep, src, fixed = TRUE)
  f <- eval(parse(text = src), envir = env); environment(f) <- env; f
}
FORM3 <- list(m = binY ~ time * trt, lv = ~time, t = ~1)
fitw <- function(d, lv) {
  warns <- character(0)
  x <- withCallingHandlers(tryCatch(mm(FORM3$m, lv.formula = lv, t.formula = FORM3$t, data = d, id = id, q = 20), error = function(e) NULL),
                           warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(x)) return(list(ok = FALSE))
  s <- summary(x)$mean.table; k <- which(rownames(s) == "time:trt")
  if (!length(k)) return(list(ok = TRUE, has_term = FALSE))
  list(ok = TRUE, has_term = TRUE, est = unname(s[k, 1]), se = unname(s[k, 2]), p = unname(s[k, 4]), code = unname(x$control["convergence_code"]),
       alpha = paste(sprintf("%.10g", x$alpha), collapse = ";"), alpha_names = paste(names(x$alpha), collapse = ";"),
       status = if (length(warns)) paste(warns, collapse = " | ") else NA_character_)
}
slope_min <- function(a, nm) { v <- as.numeric(strsplit(a, ";")[[1]]); n <- strsplit(nm, ";")[[1]]
  k <- grepl("^log\\(sigma\\):", n) & !grepl("\\(Intercept\\)", n); if (!any(k)) NA_real_ else min(exp(v[k])) }
cc_fits <- function(d, i) {
  A <- fitw(d, FORM3$lv)
  if (!A$ok) return(data.table(sim_id = i, A_ok = FALSE))
  if (!A$has_term) return(data.table(sim_id = i, A_ok = TRUE, A_has_term = FALSE))
  sm <- slope_min(A$alpha, A$alpha_names); bnd <- !(A$code %in% c(1L, 2L)) || (!is.na(sm) && sm < 1e-3)
  B <- if (bnd) fitw(d, ~1) else list(ok = NA); okB <- isTRUE(B$ok) && isTRUE(B$has_term)
  useB <- bnd && okB && is.finite(B$est) && is.finite(B$se)
  data.table(sim_id = i, A_ok = TRUE, A_has_term = TRUE, A_est = A$est, A_se = A$se, A_p = A$p, A_code = A$code, A_alpha = A$alpha,
             A_alpha_names = A$alpha_names, A_status = A$status, A_slope_min = sm, bnd = bnd, B_ok = B$ok,
             B_est = if (okB) B$est else NA_real_, B_se = if (okB) B$se else NA_real_, B_p = if (okB) B$p else NA_real_,
             B_code = if (okB) B$code else NA_integer_, B_status = if (okB) B$status else NA_character_, used = if (useB) "B" else "A",
             cc_est = if (useB) B$est else A$est, cc_se = if (useB) B$se else A$se, cc_p = if (useB) B$p else A$p,
             cc_status = if (useB) B$status else A$status)
}
# one replicate: k = seed index; gen = generator function (N, ntime, sim_id, method); analyses with model = 3
one_rep_ns <- function(env, gen, k, N, ntime, label, i) {
  deg <- design_degenerate(k, N, ntime, i)
  set_rep_seed(k, i)
  d <- tryCatch(as.data.frame(gen(N, ntime, i, label)), error = function(e) e)
  if (inherits(d, "error"))
    return(list(rows = trace_row("generation", N, ntime, label, 3L, i, paste("generation error:", conditionMessage(d))),
                cc = data.table(sim_id = i, A_ok = FALSE), rec = data.table(sim_id = i, generated = FALSE, design_degenerate = deg)))
  rows <- lapply(TECHS_ALL, function(tc) tryCatch(env$do_analysis(d, N, ntime, label, 3L, tc, i),
    error = function(e) trace_row(tc, N, ntime, label, 3L, i, paste("do_analysis error:", conditionMessage(e)))))
  list(rows = dplyr::bind_rows(rows), cc = cc_fits(d, i),
       rec = data.table(sim_id = i, generated = TRUE, design_degenerate = deg, n_obs = nrow(d), sum_y = sum(d$binY)))
}
ns_tasks <- function() {
  x <- CJ(s = seq_len(nrow(NS_SCEN)), N = as.integer(N_LEVELS), ntime = as.integer(T_LEVELS), label = c("mmlb", "manual"), chunk = 1:10)
  x <- cbind(NS_SCEN[x$s], x[, -"s"]); x[, scen := paste0("M", k)]
  x[, file := file.path(OUT_NS, "rows", sprintf("%s_m%d_N%d_T%d_%s_c%02d.rds", dgp_tag(dgp), k, N, ntime, label, chunk))]
  x
}
run_chunk_ns <- function(t, k_seed = t$k, gen_override = NULL, ids = (t$chunk - 1L) * 100L + seq_len(100L), out_file = t$file) {
  env <- load_env(t$dgp, RUNDIR_V4)
  gen <- if (!is.null(gen_override)) gen_override(env) else make_gen(env, t$dgp, t$slope)
  out <- lapply(ids, function(i) one_rep_ns(env, gen, k_seed, t$N, t$ntime, t$label, i))
  rows <- as.data.table(dplyr::bind_rows(lapply(out, `[[`, "rows")))
  rows <- rows[is.na(term) | term == "time:trt"]
  recs <- rbindlist(lapply(out, `[[`, "rec"), fill = TRUE); cc <- rbindlist(lapply(out, `[[`, "cc"), fill = TRUE)
  rows[, `:=`(dgp = t$dgp, scen = t$scen, gen_version = "v4-S6b", design_degenerate = recs$design_degenerate[match(sim_id, recs$sim_id)],
              re_structure = fifelse(method == "glmer_lrt", env$lrt_re_structure(3L), NA_character_))]
  cc[, `:=`(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
  res <- list(rows = rows, cc = cc, recs = recs)
  if (!is.null(out_file)) { dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(out_file, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out_file) }
  invisible(res)
}
