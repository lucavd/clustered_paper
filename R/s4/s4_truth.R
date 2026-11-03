# =============================================================================
# s4_truth.R -- S4 true values (D-12, D-18).
#  (1) GLMM-AR1 conditional truth per replicate: replay the design under the
#      replicate seed and solve the generator's calibration (no data generated).
#      Generator arguments (formulas, target, D, sigma2_ar, design type) are
#      CAPTURED from the v4 generator itself, not typed here.
#  (2) M1c GLMM-native truths: conditional = generator beta_c; marginal = exact
#      Gauss-Hermite quadrature on the truncated design (method of s2_c_truth.R).
#  (3) Truth table per cell for every mechanism.
# Usage: Rscript --vanilla R/s4/s4_truth.R [pilot|full]
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
source(path.expand("~/MMM/R/s3/s3_common.R"))   # load_env, set_rep_seed, gen_data, RUNDIR_V4, N_LEVELS
suppressMessages(library(data.table))
MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "full"
env <- load_env("GLMM-AR1", RUNDIR_V4)
CAL <- env$calibrate_beta_c_ar1
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
is_deg <- function(model, N, ntime, i) nrow(DEG[.(model, N, ntime, i), on = c("model", "N", "ntime", "sim_id"), nomatch = 0L]) > 0

# ---- capture the generator's calibration call ---------------------------------
gen_capture <- function(scen, model, N, ntime, label, i) {
  rec <- NULL
  env$calibrate_beta_c_ar1 <- function(...) { a <- list(...); bc <- CAL(...); rec <<- list(args = a, bc = bc); bc }
  on.exit(env$calibrate_beta_c_ar1 <- CAL)
  set_rep_seed(model, i)
  d <- tryCatch(gen_data(env, scen, model, N, ntime, i, label), error = function(e) e)
  rec
}
design_type <- function(d) if ("dose" %in% names(d)) "full" else if ("age" %in% names(d)) "age" else "none"
SC <- rbind(CJ(scen = sprintf("M%d", 1:5), label = c("mmlb", "manual"), ntime = T_LEVELS),
            CJ(scen = "M1c", label = c("mmlb", "manual"), ntime = T_M1C))
SC[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
SPEC <- lapply(seq_len(nrow(SC)), function(k) {
  s <- SC[k]; i <- 1L; while (is_deg(s$model, 25L, s$ntime, i)) i <- i + 1L
  r <- gen_capture(s$scen, s$model, 25L, s$ntime, s$label, i)
  list(f = r$args$mean.formula, re = r$args$re.formula, target = r$args$target_beta_m, D = r$args$D,
       s2 = r$args$sigma2_ar, ex = design_type(r$args$data))
})
names(SPEC) <- SC[, paste(scen, label, ntime)]
SPEC_TAB <- SC[, .(scen, label, ntime, model,
  mean_formula = vapply(SPEC, function(s) f_txt(s$f), ""), re_formula = vapply(SPEC, function(s) f_txt(s$re), ""),
  target = vapply(SPEC, function(s) num_txt(s$target), ""), D = vapply(SPEC, function(s) num_txt(s$D), ""),
  sigma2_ar = vapply(SPEC, function(s) s$s2, 0), design = vapply(SPEC, function(s) s$ex, ""))]
fwrite(SPEC_TAB, file.path(TABDIR, "V4_S4_ar1_generator_spec.csv"))

replay <- function(sp, model, N, ntime, i) {
  set_rep_seed(model, i)
  d <- env$.ar1_design(N, ntime, sp$ex)
  bc <- tryCatch(CAL(mean.formula = sp$f, re.formula = sp$re, target_beta_m = sp$target, D = sp$D,
                     sigma2_ar = sp$s2, data = d), error = function(e) conditionMessage(e))
  if (is.character(bc)) return(list(val = NA_real_, maxerr = NA_real_, it = NA_integer_, err = substr(bc, 1, 120)))
  list(val = tail(as.numeric(bc), 1), maxerr = attr(bc, "max_error"), it = attr(bc, "iterations"), err = NA_character_)
}
CELLS <- rbind(CJ(scen = sprintf("M%d", 1:5), label = c("mmlb", "manual"), N = N_LEVELS, ntime = T_LEVELS),
               CJ(scen = "M1c", label = c("mmlb", "manual"), N = N_M1C, ntime = T_M1C))
CELLS[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
reps <- 1:1000
if (MODE == "pilot") { CELLS <- CELLS[(scen == "M5" & N == 145 & ntime == 10) | (scen == "M1c" & N == 45 & ntime == 6)]; reps <- 1:20 }
t0 <- Sys.time()
one_cell <- function(k) {
  s <- CELLS[k]; sp <- SPEC[[paste(s$scen, s$label, s$ntime)]]
  r <- lapply(reps, function(i) replay(sp, s$model, s$N, s$ntime, i))
  data.table(dgp = "GLMM-AR1", scen = s$scen, N = s$N, ntime = s$ntime, label = s$label, sim_id = reps,
             truth_cond = vapply(r, `[[`, 0, "val"), max_error = vapply(r, `[[`, 0, "maxerr"),
             iterations = vapply(r, function(z) as.integer(z$it), 0L), err = vapply(r, `[[`, "", "err"))
}
TR <- rbindlist(parallel::mclapply(seq_len(nrow(CELLS)), one_cell, mc.cores = 45, mc.preschedule = FALSE))
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("AR1 replay: %d replicates in %d cells, %.1f s wall (%.4f s per replicate per core-equivalent)\n",
            nrow(TR), nrow(CELLS), el, el * min(45, nrow(CELLS)) / nrow(TR)))
if (MODE == "pilot") { print(TR[, .(n = .N, ok = sum(!is.na(truth_cond)), maxerr = max(max_error, na.rm = TRUE)), by = .(scen, label)]); quit(status = 0) }
saveRDS(TR, file.path(OUT_S4, "truth_ar1_rep.rds"), compress = FALSE)

# ---- Q6(a): reproduction of S2 in M1-M4; failures on degenerate keys only -------
S2R <- fread(file.path(REPO, "results_v4/S2/ar1_calibrated_by_replicate.csv"))[model %in% 1:4]
S2R[, scen := paste0("M", model)]
M <- merge(S2R[, .(scen, label, ntime, N, sim_id = i, s2 = beta_c)], TR[, .(scen, label, ntime, N, sim_id, v4 = truth_cond)],
           by = c("scen", "label", "ntime", "N", "sim_id"))
q6a <- M[, .(keys = .N, both_na = sum(is.na(s2) & is.na(v4)), na_mismatch = sum(is.na(s2) != is.na(v4)),
             max_abs_diff = max(abs(s2 - v4), na.rm = TRUE))]
TR[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
TR[, deg := FALSE]; TR[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
q6bij <- TR[, .(keys = .N, fail = sum(is.na(truth_cond)), deg = sum(deg), fail_and_deg = sum(is.na(truth_cond) & deg),
                fail_not_deg = sum(is.na(truth_cond) & !deg), deg_not_fail = sum(deg & !is.na(truth_cond)),
                max_error_nondeg = max(max_error[!deg], na.rm = TRUE)), by = .(block = ifelse(scen == "M1c", "M1c", "main"))]
# ---- Q6(b): replay = value used by the generator, 200 random non-degenerate keys -
set.seed(20261017L)
K <- TR[deg == FALSE][sample(.N, 200L)]
K[, gen := mapply(function(sc, m, n, t, lab, i) { r <- gen_capture(sc, m, n, t, lab, i); tail(as.numeric(r$bc), 1) },
                  scen, model, N, ntime, label, sim_id)]
q6b <- K[, .(keys = .N, max_abs_diff = max(abs(gen - truth_cond)), pass = max(abs(gen - truth_cond)) <= 1e-12)]

# ---- (2) GLMM-native truths, M1c, and the instrument check on M1 ----------------
envN <- load_env("GLMM-native", RUNDIR_V4)
gh_std <- function(q) { g <- fastGHQuad::gaussHermiteData(q); list(x = sqrt(2) * g$x, w = g$w / sqrt(pi)) }
marg_prob <- function(eta, v, gh) vapply(seq_along(eta), function(k) sum(gh$w * plogis(eta[k] + sqrt(v[k]) * gh$x)), 0)
irls <- function(X, p, w) { beta <- rep(0, ncol(X))
  for (it in 1:500) { mu <- plogis(drop(X %*% beta)); W <- w * mu * (1 - mu)
    step <- solve(crossprod(X, X * W), crossprod(X, w * (p - mu))); beta <- beta + drop(step); if (max(abs(step)) < 1e-14) break }
  beta }
marg_m1 <- function(bc, Dv, T, q = 80) {     # random intercept only, ~ time * trt, weighted truncated design
  gh <- gh_std(q); g <- CJ(time = 0:(T - 1), trt = 0:1)
  g[, w := fifelse(time == 0, 1, 0.7 + 0.3 * (T - 1 - time) / (T - 1)) * 0.5]
  X <- cbind(1, g$time, g$trt, g$time * g$trt)
  tail(irls(X, marg_prob(drop(X %*% bc), rep(Dv, nrow(X)), gh), g$w), 1)
}
PAR <- fread(file.path(TABDIR, "V4_dgm_parameters.csv"))
NT <- fread(file.path(TABDIR, "V4_S2_truth_glmmnative.csv"))
chk <- NT[model == 1, .(label, ntime, s2 = truth_marginal_weighted)]
chk[, mine := mapply(function(lab, T) marg_m1(txt_num(PAR[dgp == "GLMM-native" & model == 1 & label == lab & parameter == "beta_c"]$value),
                                             txt_mat(PAR[dgp == "GLMM-native" & model == 1 & label == lab & parameter == "D"]$value)[1, 1], T), label, ntime)]
q6c <- chk[, .(rows = .N, max_abs_diff = max(abs(s2 - mine)), pass = max(abs(s2 - mine)) <= 1e-12)]
M1C <- CJ(label = c("mmlb", "manual"), ntime = T_M1C)
M1C[, truth_conditional := mapply(function(lab, T) tail(envN$M1C_ZEGER_F * envN$m1c_beta_marginal(T, lab), 1), label, ntime)]
M1C[, truth_marginal := mapply(function(lab, T) marg_m1(envN$M1C_ZEGER_F * envN$m1c_beta_marginal(T, lab), envN$M1C_VAR, T), label, ntime)]
M1C[, quad_err_40v80 := mapply(function(lab, T) abs(marg_m1(envN$M1C_ZEGER_F * envN$m1c_beta_marginal(T, lab), envN$M1C_VAR, T, 40) -
                                                    marg_m1(envN$M1C_ZEGER_F * envN$m1c_beta_marginal(T, lab), envN$M1C_VAR, T, 80)), label, ntime)]
envG <- load_env("GenBinaryY", RUNDIR_V4)

# ---- (3) truth table per (dgp, scen, ntime, label) -----------------------------
TV <- fread(file.path(TABDIR, "V4_S2_truth_values.csv"))
gb <- TV[dgp == "GenBinaryY", .(dgp, scen = paste0("M", model), label, truth_marginal)]
gb <- gb[CJ(scen = sprintf("M%d", 1:5), ntime = T_LEVELS, label = c("mmlb", "manual")), on = c("scen", "label")]
gb <- rbind(gb, CJ(dgp = "GenBinaryY", scen = "M1c", ntime = T_M1C, label = c("mmlb", "manual"))[,
  truth_marginal := mapply(function(lab, T) tail(envG$m1c_beta_marginal(T, lab), 1), label, ntime)], use.names = TRUE)
gb[, `:=`(truth_conditional = NA_real_, cond_source = "none (D-18)", marg_source = "generator beta")]
nt <- rbind(NT[, .(dgp, scen = paste0("M", model), ntime, label, truth_marginal = truth_marginal_weighted, truth_conditional,
                   cond_source = "generator beta_c", marg_source = "GH quadrature, weighted design (S2)")],
            M1C[, .(dgp = "GLMM-native", scen = "M1c", ntime, label, truth_marginal, truth_conditional,
                    cond_source = "generator beta_c", marg_source = "GH quadrature, weighted design (S4)")])
ar <- SC[, .(dgp = "GLMM-AR1", scen, ntime, label,
             truth_marginal = vapply(SPEC, function(s) tail(as.numeric(s$target), 1), 0), truth_conditional = NA_real_,
             cond_source = "per replicate (results_v4/S4/truth_ar1_rep.rds)", marg_source = "generator target")]
TT <- rbindlist(list(gb, nt, ar), use.names = TRUE)
setcolorder(TT, c("dgp", "scen", "ntime", "label"))
fwrite(TT[order(dgp, scen, ntime, label)], file.path(TABDIR, "V4_S4_truth_by_cell.csv"))
fwrite(M1C, file.path(TABDIR, "V4_S4_truth_m1c_native.csv"))
Q6 <- list(a = q6a, bij = q6bij, b = q6b, c = q6c, spec = SPEC_TAB)
saveRDS(Q6, file.path(TABDIR, "V4_S4_q6.rds"))
fwrite(rbindlist(list(
  data.table(check = "Q6a S2 reproduction M1-M4", keys = q6a$keys, value = q6a$max_abs_diff, na_mismatch = q6a$na_mismatch,
             pass = q6a$max_abs_diff <= 1e-12 & q6a$na_mismatch == 0),
  q6bij[, .(check = paste("Q6a failures = degenerate keys,", block), keys, value = max_error_nondeg, na_mismatch = fail_not_deg + deg_not_fail,
            pass = fail_not_deg == 0 & deg_not_fail == 0 & max_error_nondeg < 1e-8)],
  data.table(check = "Q6b replay = generator, 200 keys", keys = q6b$keys, value = q6b$max_abs_diff, na_mismatch = 0L, pass = q6b$pass),
  data.table(check = "Q6c quadrature reproduces S2 M1", keys = q6c$rows, value = q6c$max_abs_diff, na_mismatch = 0L, pass = q6c$pass))),
  file.path(TABDIR, "V4_S4_q6_truth_checks.csv"))
cat("Q6 checks:\n"); print(fread(file.path(TABDIR, "V4_S4_q6_truth_checks.csv")))
cat("M1c GLMM-native truths:\n"); print(M1C, digits = 6)
cat("AR1 conditional truth per cell (mean, sd) at N = 145 / M1c N = 45:\n")
print(dcast(TR[deg == FALSE & (N == 145 | (scen == "M1c" & N == 45)), .(m = sprintf("%.4f (%.1e)", mean(truth_cond), sd(truth_cond))), by = .(scen, label, ntime)],
            scen + label ~ ntime, value.var = "m"))
