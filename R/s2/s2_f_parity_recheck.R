# =============================================================================
# s2_f_parity_recheck.R -- S2(f), Q6.1 re-check. The first check (in
# s2_f_pilot.R) compared GenBinaryY arguments with identical(), which for a
# formula also compares its environment; it flagged mean/lv/t.formula as
# different even between the ORIGINAL H0 and H1, whose source text is the same.
# Here formulas are compared by deparse(), data by value (data.frame columns),
# numerics exactly. Generators: .rundir (original) and R/s2/pilot_gen (v4 H0).
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
GEN_V4 <- file.path(REPO, "R", "s2", "pilot_gen")
CELLS <- CJ(model = 3:5, NT = 1:4)[, `:=`(N = c(10, 25, 65, 145)[NT], ntime = c(2, 4, 6, 10)[NT])][, NT := NULL]
mk_env <- function(v4) { env <- load_family("GenBinaryY", RUNDIR)
  if (v4) for (f in sprintf("generation_model%d.R", 3:5)) sys.source(file.path(GEN_V4, f), envir = env); env }
E_OLD <- mk_env(FALSE); E_V4 <- mk_env(TRUE)
capture_args <- function(env, m, N, Tt, lab, i) {
  got <- NULL; real <- MMLB::GenBinaryY
  env$GenBinaryY <- function(...) { got <<- list(...); real(...) }
  on.exit(rm("GenBinaryY", envir = env)); set_rep_seed(m, i); env$make_all_data(N, Tt, i, m, lab); got
}
canon <- function(x) {
  if (inherits(x, "formula")) paste(deparse(x, width.cutoff = 500L), collapse = "")
  else if (is.data.frame(x)) lapply(as.list(x), as.vector) else x
}
same <- function(a, b) identical(canon(a), canon(b))
diffs <- function(a, b) paste(names(a)[!vapply(names(a), function(n) same(a[[n]], b[[n]]), TRUE)], collapse = ",")
par <- rbindlist(lapply(seq_len(nrow(CELLS)), function(k) {
  cl <- CELLS[k]
  a1 <- capture_args(E_OLD, cl$model, cl$N, cl$ntime, "mmlb", 1)
  a0 <- capture_args(E_V4, cl$model, cl$N, cl$ntime, "manual", 1)
  a0o <- capture_args(E_OLD, cl$model, cl$N, cl$ntime, "manual", 1)
  b_pos <- paste(which(a1$beta != a0$beta), collapse = ",")
  data.table(cl, v4H0_vs_H1_differs_in = diffs(a1, a0), beta_diff_positions = b_pos, beta_len = length(a1$beta),
             origH0_vs_H1_differs_in = diffs(a1, a0o),
             sigma_H1 = paste(a1$sigma, collapse = ","), sigma_v4H0 = paste(a0$sigma, collapse = ","),
             gamma_H1 = paste(a1$gamma, collapse = ","), gamma_v4H0 = paste(a0$gamma, collapse = ","))
}))
par[, pass := v4H0_vs_H1_differs_in == "beta" & beta_diff_positions == as.character(beta_len)]
options(width = 250); print(par)
cat(sprintf("Q6.1 parity (re-check, formulas by deparse): %d/%d cells\n", sum(par$pass), nrow(par)))
# controprova of the instrument: original H0 must differ exactly as the S1 table says (P2)
exp_orig <- c("3" = "beta,sigma", "4" = "beta,sigma,gamma", "5" = "beta,sigma,gamma")
par[, instr_ok := origH0_vs_H1_differs_in == exp_orig[as.character(model)]]
cat(sprintf("Q6.1 instrument controprova (original H0 differs as in S1 P2): %d/%d cells\n", sum(par$instr_ok), nrow(par)))
fwrite(par, file.path(TABDIR, "V4_S2_pilot_parity.csv"))
