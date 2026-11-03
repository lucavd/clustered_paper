
# s10_lead_sep_check.R -- S10, lead check of the Degeneri result at the two flagged points of GEE AR(1) Fay-Graubard
# (M2, N = 10, T = 2, H0, GLMM-AR1 and GLMM-native), from the raw S7d rows and an LP separation check written here.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L)
source(path.expand("~/MMM/R/s7d/s7d_common.R")); .libPaths(c(.libPaths(), path.expand("~/R/s10_lib"))); suppressMessages(library(lpSolveAPI))
CL <- readRDS("/home/user/MMM/results_v4/S10/degen/classes.rds"); cat("classes columns:", paste(names(CL), collapse = ","), "\n")
K <- as.data.table(readRDS("/home/user/MMM/results_v4/S9/const_keys.rds"))
sep_lp <- function(X, y) { A <- X * (2 * y - 1); p <- ncol(A); lp <- make.lp(nrow(A), p)
  for (j in seq_len(p)) set.column(lp, j, A[, j]); set.constr.type(lp, rep(">=", nrow(A))); set.rhs(lp, rep(0, nrow(A)))
  set.bounds(lp, lower = rep(-1, p), upper = rep(1, p)); set.objfn(lp, colSums(A)); lp.control(lp, sense = "max"); solve(lp); get.objective(lp) > 1e-7 }
TK <- s7d_tasks()
for (g in c("GLMM-AR1", "GLMM-native")) {
  tk <- TK[dgp == g & scen == "M2" & N == 10 & ntime == 2 & label == "manual"]
  env <- load_env(g, RUNDIR_V4); gen <- s7d_gen(env, tk[1])
  out <- rbindlist(lapply(seq_len(nrow(tk)), function(r) { x <- readRDS(tk$file[r]); rec <- as.data.table(x$recs); rw <- as.data.table(x$rows)[cor == "ar1"]
    merge(rw[, .(sim_id, gee_ok, est, se_fg)], rec[, .(sim_id, design_degenerate, n_obs, sum_y)], by = "sim_id") }))
  kk <- K[dgp == g & scen == "M2" & N == 10 & ntime == 2 & label == "manual", sim_id]
  out[, use := !design_degenerate & gee_ok & is.finite(est) & abs(est) < 10 & is.finite(se_fg) & se_fg > 0 & se_fg < 10 & !(sim_id %in% kk)]
  out[, rej := use & 2 * pt(-abs(est / se_fg), 10 - 5L) < 0.05]
  out[, sep := NA]
  for (i in out[use == TRUE, sim_id]) { set_rep_seed(tk$model[1], i); d <- as.data.frame(gen(10L, 2L, i, "manual"))
    out[sim_id == i, sep := sep_lp(model.matrix(FORM_OF("M2"), d), d$binY)] }
  cat(sprintf("%s: usable %d, rejections %d | usable separated %d, rejections on separated %d | usable not separated %d, rejections %d\n", g,
      out[use == TRUE, .N], out[rej == TRUE, .N], out[use & sep %in% TRUE, .N], out[rej & sep %in% TRUE, .N], out[use & sep %in% FALSE, .N], out[rej & sep %in% FALSE, .N]))
}
