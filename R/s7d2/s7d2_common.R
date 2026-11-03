# =============================================================================
# s7d2_common.R -- S7d2: interval by inversion of the sign-flip score test on the 2,376 design points of the four
# mechanisms (D-56): the 1,896 points of the S7d run (R/s7d/s7d_common.R::s7d_tasks) and the 480 copula points of the
# S7e run (R/s7e/s7e_common.R::s7e_tasks); H0 and H1; 1,000 replicates per arm; 47,520 chunks of 100 replicates.
# Datasets regenerated exactly as in those runs (degenerate-design flag, set_rep_seed(model, i), same generator); then
#   set.seed(20261022 + i); s7d_flip()            -> p0_old (reproduction of the saved p_flip, gate GG2)
#   set.seed(20261022 + i); flip_ci(d, scen, th)  -> R/s7d2/s7d2_ci.R (p at 0 and at the true value, centre, interval)
# The saved S7d / S7e chunk of the same replicates is read to carry n_obs, sum_y, p_flip (GG2) and the GEE AR(1)
# with fewclust Fay-Graubard (est, se_fg, gee_ok) used as comparator.
# =============================================================================
source(path.expand("~/MMM/R/s7e/s7e_common.R"))   # s7d_common (s3_common, s6b_ns_common, s7c_tests, s7d_flip), copula env
source(path.expand("~/MMM/R/s7d2/s7d2_ci.R"))
OUT_S7D2 <- file.path(REPO, "results_v4/S7d2"); dir.create(file.path(OUT_S7D2, "rows"), recursive = TRUE, showWarnings = FALSE)
TRUTH <- rbind(fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))[, .(dgp, scen, ntime, label, th = truth_marginal)],
               fread(file.path(TABDIR, "V4_S6b_ns_truth.csv"))[, .(dgp, scen, ntime, label, th = tm)])
truth_of <- function(g, s, nt, lb) {
  if (g == DGP_COP) return(if (lb == "manual") 0 else if (s == "M4") 0.05 else 0.2)
  x <- TRUTH[dgp == g & scen == s & ntime == nt & label == lb]
  stopifnot(nrow(x) == 1L, is.finite(x$th)); x$th
}
s7d2_tasks <- function() {
  a <- s7d_tasks(); a[, old_file := file]
  b <- s7e_tasks(); b[, `:=`(dgp = DGP_COP, old_file = file)]
  x <- rbind(a[, .(dgp, scen, N, ntime, label, chunk, model, old_file)], b[, .(dgp, scen, N, ntime, label, chunk, model, old_file)])
  stopifnot(nrow(x) == 47520L, nrow(unique(x[, .(dgp, scen, N, ntime)])) == 2376L)
  x[, file := file.path(OUT_S7D2, "rows", basename(old_file))]
  x
}
s7d2_gen <- function(t) {
  if (t$dgp == DGP_COP) { env <- load_env_cop(); function(N, ntime, i, label) env$make_all_data(N, ntime, i, t$model, label) }
  else { env <- load_env(t$dgp, RUNDIR_V4); s7d_gen(env, t) }
}
CI_COLS <- c("p0_new", "p0_guard", "why0_guard", "p_theta", "theta_why", "est", "se_mod", "est_why", "lo", "hi", "rej_lo", "rej_hi", "p_lo", "p_hi",
             "p_lo_out", "p_hi_out", "p_est", "su_est", "n_eval", "n_na")
one_rep_s7d2 <- function(gen, t, i, th) {
  deg <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) e)
  rec <- data.table(sim_id = i, generated = !inherits(d, "error"), design_degenerate = deg, n_obs = NA_integer_, sum_y = NA_real_,
                    theta = th, p0_old = NA_real_, ci_error = NA_character_, secs = NA_real_)
  if (inherits(d, "error")) return(rec)
  d <- d[order(d$id, d$time), ]; rec$n_obs <- nrow(d); rec$sum_y <- sum(d$binY)
  set.seed(FLIP_SEED(i)); fl <- tryCatch(s7d_flip(d, t$scen), error = function(e) NULL)
  if (!is.null(fl)) rec$p0_old <- fl$p[1]
  set.seed(FLIP_SEED(i)); t0 <- proc.time()[[3]]
  ci <- tryCatch(flip_ci(d, t$scen, th), error = function(e) e)
  rec$secs <- proc.time()[[3]] - t0
  if (inherits(ci, "error")) { rec$ci_error <- substr(conditionMessage(ci), 1, 160); return(rec) }
  cbind(rec, as.data.table(ci[CI_COLS]))
}
old_values <- function(t) {
  o <- readRDS(t$old_file)
  g <- if (t$dgp == DGP_COP) o$gee else o$rows
  merge(o$recs[, .(sim_id, old_n_obs = n_obs, old_sum_y = sum_y, old_p_flip = p_flip)],
        g[cor == "ar1", .(sim_id, ar1_gee_ok = gee_ok, ar1_est = est, ar1_se_fg = se_fg)], by = "sim_id", all.x = TRUE)
}
run_chunk_s7d2 <- function(t, ids = (t$chunk - 1L) * 100L + seq_len(100L), out_file = t$file) {
  gen <- s7d2_gen(t); th <- truth_of(t$dgp, t$scen, t$ntime, t$label)
  recs <- rbindlist(lapply(ids, function(i) one_rep_s7d2(gen, t, i, th)), fill = TRUE)
  recs <- merge(recs, old_values(t), by = "sim_id", all.x = TRUE)
  recs[, `:=`(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label)]
  if (!is.null(out_file)) { tmp <- paste0(out_file, ".tmp"); saveRDS(recs, tmp); file.rename(tmp, out_file) }
  invisible(recs)
}
