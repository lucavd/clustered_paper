.libPaths("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
Sys.setenv(OPENBLAS_NUM_THREADS="1")
suppressMessages({library(MMLB); library(lme4); library(geepack); library(geesmv); library(dplyr); library(tibble); library(MASS); library(fastGHQuad)})
source("R/protocol/dgm_common.R"); RNGkind("L'Ecuyer-CMRG","Inversion","Rejection")
env <- load_family("GLMM-native", ".rundir")
for (f in c(sprintf("analysis_model%d.R",1:5),"analysis_smallsample.R","do_analysis.R")) sys.source(file.path(".rundir",f), envir=env)
for (m in c(3,5)) for (cell in list(c(145,10), c(25,4))) {
  N <- cell[1]; Tt <- cell[2]; tt <- c(ri=0, rs=0); nsing <- c(ri=0, rs=0)
  re <- if (m == 5) "(1 + time + dose | id)" else "(1 + time | id)"
  rhs <- if (m == 5) "time*trt + dose" else "time*trt"
  for (i in 1:5) {
    set.seed(seed_of(m, i)); d <- env$make_all_data(N, Tt, i, m, "mmlb")
    for (k in c("ri","rs")) {
      f <- as.formula(paste("binY ~", rhs, "+", if (k=="ri") "(1|id)" else re))
      s <- system.time(r <- suppressWarnings(suppressMessages(env$analyze_glmer_lrt(d, f, N, Tt, "mmlb"))))[["elapsed"]]
      tt[k] <- tt[k] + s; nsing[k] <- nsing[k] + is.null(r)
    }
  }
  cat(sprintf("M%d N=%d T=%d: s/replicate LRT (1|id) %.2f, slopes %.2f | NULL results ri %d/5, rs %d/5\n", m, N, Tt, tt["ri"]/5, tt["rs"]/5, nsing["ri"], nsing["rs"]))
}
