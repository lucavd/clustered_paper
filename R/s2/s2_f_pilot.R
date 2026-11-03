# =============================================================================
# s2_f_pilot.R -- S2(f), Q6: pilot of the corrected H0 arm (D-01) for
# GenBinaryY M3-M5 on 12 cells x 1000 replicates. Same seed, explicit RNGkind
# (D-07), all techniques of the executed code. Output: results_v4/pilot_S2/.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
PILOT <- file.path(REPO, "results_v4", "pilot_S2"); dir.create(PILOT, recursive = TRUE, showWarnings = FALSE)
GEN_V4 <- file.path(REPO, "R", "s2", "pilot_gen"); dir.create(GEN_V4, showWarnings = FALSE)
TECHS <- c("mmm", "gee", "glmer", "gee_ss", "glmer_lrt")
CELLS <- CJ(model = 3:5, NT = 1:4)[, `:=`(N = c(10, 25, 65, 145)[NT], ntime = c(2, 4, 6, 10)[NT])][, NT := NULL]
ANALYSIS_FILES <- c(sprintf("analysis_model%d.R", 1:5), "analysis_smallsample.R", "do_analysis.R")

# ---- corrected H0 generators: only the H0 sigma / gamma lines change ----------
EDITS <- list(
  "generation_model3.R" = list(c(34, "    manual = c(1, 0))", "    manual = c(1, 0.01))")),
  "generation_model4.R" = list(c(35, "    manual = c(1, 0) )", "    manual = c(1, 0.01) ) # D-01: H0 sigma = H1 sigma"),
                               c(40, "    manual = c(1, 0) )", "    manual = c(1, 0.01) )")),
  "generation_model5.R" = list(c(42, "    manual = c(1, 0, 0))", "    manual = c(1, 0.01, 0.01))"),
                               c(47, "     manual = c(1, 0, 0))", "     manual = c(1, 0.1, 0.05))")))
for (f in names(EDITS)) {
  L <- readLines(file.path(RUNDIR, f), encoding = "UTF-8")
  for (e in EDITS[[f]]) { k <- as.integer(e[1]); stopifnot(startsWith(L[k], e[2])); L[k] <- e[3] }
  writeLines(c(sprintf("# PILOT S2 (D-01): copy of .rundir/%s with H0 sigma/gamma set to the H1 values.", f), L),
             file.path(GEN_V4, f), useBytes = TRUE)
}
mk_env <- function(v4) {
  env <- load_family("GenBinaryY", RUNDIR)
  if (v4) for (f in names(EDITS)) sys.source(file.path(GEN_V4, f), envir = env)
  for (f in ANALYSIS_FILES) sys.source(file.path(RUNDIR, f), envir = env)
  env
}
E_OLD <- mk_env(FALSE); E_V4 <- mk_env(TRUE)

# ---- Q6.1 parity: arguments handed to GenBinaryY ---------------------------------
capture_args <- function(env, m, N, Tt, lab, i) {
  got <- NULL; real <- MMLB::GenBinaryY
  env$GenBinaryY <- function(...) { got <<- list(...); real(...) }
  on.exit(rm("GenBinaryY", envir = env))
  set_rep_seed(m, i); env$make_all_data(N, Tt, i, m, lab); got
}
par <- rbindlist(lapply(seq_len(nrow(CELLS)), function(k) {
  cl <- CELLS[k]; a1 <- capture_args(E_OLD, cl$model, cl$N, cl$ntime, "mmlb", 1)
  a0 <- capture_args(E_V4, cl$model, cl$N, cl$ntime, "manual", 1)
  a0o <- capture_args(E_OLD, cl$model, cl$N, cl$ntime, "manual", 1)
  others <- setdiff(names(a1), "beta")
  data.table(cl, same_names = identical(names(a1), names(a0)),
             others_identical = all(vapply(others, function(n) identical(a1[[n]], a0[[n]]), TRUE)),
             beta_diff_positions = paste(which(a1$beta != a0$beta), collapse = ","), beta_len = length(a1$beta),
             old_H0_differs_in = paste(setdiff(others[!vapply(others, function(n) identical(a1[[n]], a0o[[n]]), TRUE)], "data"), collapse = ","))
}))
par[, pass := same_names & others_identical & beta_diff_positions == as.character(beta_len)]
print(par); cat(sprintf("Q6.1 parity: %d/%d cells\n", sum(par$pass), nrow(par)))
fwrite(par, file.path(TABDIR, "V4_S2_pilot_parity.csv"))

# ---- one replicate, exactly as run_chunk() ----------------------------------------
one_rep <- function(env, m, N, Tt, lab, i) {
  set_rep_seed(m, i)
  d <- tryCatch(env$make_all_data(N, Tt, i, m, lab), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  rows <- lapply(TECHS, function(tc) tryCatch(env$do_analysis(d, N, Tt, lab, m, tc, i), error = function(e) NULL))
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)) dplyr::bind_rows(rows) else NULL
}

# ---- Q6.2 controprova: the bench reproduces results_full for i = 1..5 --------------
rep_chk <- rbindlist(parallel::mcMap(function(k, i) {
  cl <- CELLS[k]
  st <- as.data.table(readRDS(sprintf("%s/results_full/GenBinaryY_m%d_N%d_T%d_manual_c01.rds", REPO, cl$model, cl$N, cl$ntime)))[sim_id == i]
  rf <- as.data.table(one_rep(E_OLD, cl$model, cl$N, cl$ntime, "manual", i))
  key <- function(x) paste(x$method, x$term, seq_along(x$term))
  same <- nrow(st) == nrow(rf) && identical(key(st), key(rf)) &&
    identical(st$estimate, rf$estimate) && identical(st$se, rf$se) && identical(st$p, rf$p) && identical(st$status, rf$status)
  data.table(cl, i = i, rows_stored = nrow(st), rows_bench = nrow(rf), identical = same,
             max_d_est = if (nrow(st) == nrow(rf)) suppressWarnings(max(abs(st$estimate - rf$estimate), na.rm = TRUE)) else NA_real_)
}, rep(seq_len(nrow(CELLS)), each = 5), rep(1:5, nrow(CELLS)), mc.cores = 40))
print(rep_chk[identical == FALSE]); cat(sprintf("Q6.2 bench reproduces results_full: %d/%d\n", sum(rep_chk$identical), nrow(rep_chk)))
fwrite(rep_chk, file.path(TABDIR, "V4_S2_pilot_reproduction.csv"))

# ---- pilot run: corrected H0, 12 cells x 10 chunks x 100 ---------------------------
tasks <- CELLS[, .(chunk = 1:10), by = .(model, N, ntime)]
tasks[, file := sprintf("%s/GenBinaryY_m%d_N%d_T%d_manualv4_c%02d.rds", PILOT, model, N, ntime, chunk)]
tasks <- tasks[order(-N * ntime)]          # longest first
t0 <- Sys.time()
res <- parallel::mcMap(function(m, N, Tt, ch, f) {
  if (file.exists(f)) return("skip")
  acc <- list(); des <- list()
  for (j in 1:100) {
    i <- (ch - 1L) * 100L + j
    r <- one_rep(E_V4, m, N, Tt, "manual", i)
    if (!is.null(r)) acc[[length(acc) + 1L]] <- r
    set_rep_seed(m, i); dn <- tryCatch(E_V4$make_all_data(N, Tt, i, m, "manual"), error = function(e) NULL)
    set_rep_seed(m, i); do <- tryCatch(E_OLD$make_all_data(N, Tt, i, m, "manual"), error = function(e) NULL)
    cols <- intersect(c("id", "time", "trt", "dose"), names(dn))
    des[[j]] <- data.table(sim_id = i, design_identical = !is.null(dn) && !is.null(do) && identical(as.data.frame(dn)[cols], as.data.frame(do)[cols]),
                           y_identical_share = if (!is.null(dn) && !is.null(do)) mean(dn$binY == do$binY) else NA_real_)
  }
  out <- dplyr::bind_rows(acc); out$dgp <- "GenBinaryY"; out$label_v4 <- "manual_v4"
  attr(out, "design_check") <- rbindlist(des)
  saveRDS(out, paste0(f, ".tmp"), compress = "xz"); file.rename(paste0(f, ".tmp"), f); "done"
}, tasks$model, tasks$N, tasks$ntime, tasks$chunk, tasks$file, mc.cores = 40, mc.preschedule = FALSE)
cat(sprintf("pilot: %s | %.1f min\n", paste(names(table(unlist(res))), table(unlist(res)), collapse = " "),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
errs <- which(!vapply(res, function(z) is.character(z) && z %in% c("done", "skip"), TRUE))
if (length(errs)) { cat("task errors:\n"); print(res[errs]) }
