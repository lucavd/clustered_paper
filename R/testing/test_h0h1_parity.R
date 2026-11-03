# =============================================================================
# test_h0h1_parity.R -- S3 Q1/Q2. For every generator (3 mechanisms x M1-M5,
# plus M1c when present) the arguments handed to the generating engine under H0
# and H1 must differ ONLY in the last element of the coefficient vector
# (GenBinaryY: beta; GLMM-native: beta_c; GLMM-AR1: the calibration target, with
# every other argument of calibrate_beta_c_ar1 and gen_ar1_data identical).
# Formulas compared by deparse(), data by value, numbers exactly.
# Usage: Rscript --vanilla R/testing/test_h0h1_parity.R [rundir]   (default R/v4/rundir)
# =============================================================================
source(path.expand("~/MMM/R/s3/s3_common.R"))
a <- commandArgs(TRUE); RD <- if (length(a)) normalizePath(a[1]) else RUNDIR_V4
canon <- function(x) {
  if (inherits(x, "formula")) paste(deparse(x, width.cutoff = 500L), collapse = "")
  else if (is.data.frame(x)) lapply(as.list(x), as.vector) else x
}
same  <- function(a, b) identical(canon(a), canon(b))
diffs <- function(a, b) names(a)[!vapply(names(a), function(n) same(a[[n]], b[[n]]), TRUE)]
capture <- function(env, dgp, scen, model, N, Tt, lab, i) {
  got <- list()
  if (dgp == "GenBinaryY") {
    real <- MMLB::GenBinaryY; env$GenBinaryY <- function(...) { got$engine <<- list(...); real(...) }
    on.exit(rm("GenBinaryY", envir = env))
  } else if (dgp == "GLMM-native") {
    real <- env$gen_glmm_data; env$gen_glmm_data <- function(...) { got$engine <<- list(...); real(...) }
    on.exit(assign("gen_glmm_data", real, envir = env))
  } else {
    rc <- env$calibrate_beta_c_ar1; rg <- env$gen_ar1_data
    env$calibrate_beta_c_ar1 <- function(...) { got$calib <<- list(...); rc(...) }
    env$gen_ar1_data <- function(...) { x <- list(...); x$beta_c <- NULL; got$engine <<- x; rg(...) }
    on.exit({ assign("calibrate_beta_c_ar1", rc, envir = env); assign("gen_ar1_data", rg, envir = env) })
  }
  # arguments are recorded before the engine runs; a generation failure (AR1
  # calibration on a degenerate design, S2) does not prevent the comparison
  set_rep_seed(model, i); tryCatch(gen_data(env, scen, model, N, Tt, i, lab), error = function(e) NULL); got
}
# gen_ar1_data is not reached when the calibration fails (degenerate design)
gen_part <- function(h) if (is.null(h$engine)) list() else setNames(h$engine, paste0("gen.", names(h$engine)))
COEF <- c("GenBinaryY" = "beta", "GLMM-native" = "beta_c", "GLMM-AR1" = "target_beta_m")
CELLS_MAIN <- list(c(10, 2), c(25, 4), c(45, 6), c(145, 10)); CELLS_M1C <- list(c(10, 2), c(25, 4), c(45, 6))
has_m1c <- file.exists(file.path(RD, "m1c_params.R"))
gens <- expand.grid(dgp = DGPS, scen = c(sprintf("M%d", 1:5), if (has_m1c) "M1c"), stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(gens)), function(g) {
  dgp <- gens$dgp[g]; scen <- gens$scen[g]; model <- if (scen == "M1c") 1L else as.integer(sub("M", "", scen))
  env <- load_env(dgp, RD); cells <- if (scen == "M1c") CELLS_M1C else CELLS_MAIN
  chk <- do.call(rbind, lapply(cells, function(cl) do.call(rbind, lapply(1:3, function(i) {
    h1 <- capture(env, dgp, scen, model, cl[1], cl[2], "mmlb", i)
    h0 <- capture(env, dgp, scen, model, cl[1], cl[2], "manual", i)
    if (dgp == "GLMM-AR1") { a1 <- c(h1$calib, gen_part(h1)); a0 <- c(h0$calib, gen_part(h0)) }
    else { a1 <- h1$engine; a0 <- h0$engine }
    dd <- diffs(a1, a0); cf <- COEF[[dgp]]
    pos <- which(a1[[cf]] != a0[[cf]])
    data.frame(N = cl[1], ntime = cl[2], i = i, differs_in = paste(dd, collapse = ","),
               coef_positions = paste(pos, collapse = ","), coef_len = length(a1[[cf]]),
               pass = identical(dd, cf) && identical(pos, length(a1[[cf]])))
  }))))
  data.frame(dgp = dgp, scen = scen, n_checks = nrow(chk), n_pass = sum(chk$pass),
             differs_in = paste(unique(chk$differs_in), collapse = " / "), pass = all(chk$pass))
}))
options(width = 200); print(res)
main <- res$scen != "M1c"
cat(sprintf("PARITY [%s]: main grid %d/%d PASS%s\n", RD, sum(res$pass[main]), sum(main),
            if (has_m1c) sprintf(" | M1c %d/%d PASS", sum(res$pass[!main]), sum(!main)) else ""))
if (nzchar(Sys.getenv("MMM_PARITY_OUT"))) write.csv(res, Sys.getenv("MMM_PARITY_OUT"), row.names = FALSE)

# ---- Q2 (D-02): AR1 M5 v4 vs original, and D equal to GLMM-native M5 ---------
if (has_m1c) {
  ev4 <- load_env("GLMM-AR1", RD); eor <- load_env("GLMM-AR1", RUNDIR_ORIG); egn <- load_env("GLMM-native", RD)
  q2 <- do.call(rbind, lapply(CELLS_MAIN, function(cl) do.call(rbind, lapply(c("mmlb", "manual"), function(lab) {
    v <- capture(ev4, "GLMM-AR1", "M5", 5L, cl[1], cl[2], lab, 1); o <- capture(eor, "GLMM-AR1", "M5", 5L, cl[1], cl[2], lab, 1)
    gn <- capture(egn, "GLMM-native", "M5", 5L, cl[1], cl[2], lab, 1)
    dc <- diffs(v$calib, o$calib); de <- diffs(v$engine, o$engine)
    Dv <- as.matrix(v$calib$D); Do <- as.matrix(o$calib$D)
    data.frame(N = cl[1], ntime = cl[2], label = lab, calib_differs = paste(dc, collapse = ","),
               gen_differs = paste(de, collapse = ","), D_cells_changed = paste(which(Dv != Do), collapse = ","),
               D33_v4 = Dv[3, 3], D33_orig = Do[3, 3], D_equal_glmmnative = identical(Dv, as.matrix(gn$engine$D)),
               D_gen_equal_calib = identical(Dv, as.matrix(v$engine$D)))
  }))))
  q2$pass <- q2$calib_differs == "D" & q2$gen_differs == "D" & q2$D_cells_changed == "9" &
    q2$D33_v4 == 0.001 & q2$D33_orig == 0.01 & q2$D_equal_glmmnative & q2$D_gen_equal_calib
  print(q2); cat(sprintf("Q2 D-02: %d/%d PASS\n", sum(q2$pass), nrow(q2)))
  if (nzchar(Sys.getenv("MMM_Q2_OUT"))) write.csv(q2, Sys.getenv("MMM_Q2_OUT"), row.names = FALSE)
}
