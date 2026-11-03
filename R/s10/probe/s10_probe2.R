
# s10_probe2.R -- S10 probe 2: direct LP separation check (Konis 2007) vs detectseparation, timings incl. largest cell
source(path.expand("~/MMM/R/s7d/s7d_common.R"))
.libPaths(c(.libPaths(), path.expand("~/R/s10_lib")))
suppressMessages({ library(detectseparation); library(lpSolveAPI) })
# MLE exists iff no beta != 0 with s_i x_i'beta >= 0 for all i (s = 2y - 1). LP: max sum(s_i x_i'beta) s.t. s_i x_i'beta >= 0, |beta_j| <= 1.
lp_sep <- function(X, y) {
  A <- X * (2 * y - 1); p <- ncol(A)
  lp <- make.lp(nrow(A), p); for (j in seq_len(p)) set.column(lp, j, A[, j])
  set.constr.type(lp, rep(">=", nrow(A))); set.rhs(lp, rep(0, nrow(A)))
  set.bounds(lp, lower = rep(-1, p), upper = rep(1, p)); set.objfn(lp, colSums(A)); lp.control(lp, sense = "max")
  st <- solve(lp); if (st != 0) return(NA)
  get.objective(lp) > 1e-7
}
TK <- s7d_tasks()
pick <- TK[(dgp == "GenBinaryY" & scen == "M2" & N == 10 & ntime == 2 & label == "manual" & chunk == 1) |
           (dgp == "GLMM-AR1" & scen == "M5" & N == 25 & ntime == 3 & label == "mmlb" & chunk == 1) |
           (dgp == "GLMM-native" & scen == "M3" & N == 15 & ntime == 4 & label == "manual" & chunk == 1) |
           (dgp == "GenBinaryY" & scen == "M4" & N == 145 & ntime == 10 & label == "mmlb" & chunk == 1)]
for (r in seq_len(nrow(pick))) { t <- pick[r]; env <- load_env(t$dgp, RUNDIR_V4); gen <- s7d_gen(env, t)
  tg <- 0; tl <- 0; td <- 0; agree <- 0L; n <- 0L; sl <- 0L
  for (i in 1:100) { set_rep_seed(t$model, i); t0 <- proc.time()[[3]]
    d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) NULL); tg <- tg + proc.time()[[3]] - t0; if (is.null(d)) next
    X <- model.matrix(FORM_OF(t$scen), d); y <- d$binY
    t0 <- proc.time()[[3]]; a <- lp_sep(X, y); tl <- tl + proc.time()[[3]] - t0
    t0 <- proc.time()[[3]]; s <- glm(FORM_OF(t$scen), data = d, family = binomial, method = "detect_separation"); td <- td + proc.time()[[3]] - t0
    n <- n + 1L; agree <- agree + identical(isTRUE(a), isTRUE(s$outcome)); sl <- sl + isTRUE(a) }
  cat(sprintf("%s %s N%d T%d %s: gen %.2fs, LP %.2fs, detectsep %.2fs per 100; separated(LP) %d; agree %d/%d\n",
              t$dgp, t$scen, t$N, t$ntime, t$label, tg, tl, td, sl, agree, n)) }
