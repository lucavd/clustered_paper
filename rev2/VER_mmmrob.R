.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
# VER_mmmrob.R -- verifica preliminare: mean.table di summary.MMLongit e' mod.cov,
# rob.cov e' presente e diversa. Poi calibrazione SE/SD a N=60, T=6 (M-30).
suppressMessages({library(data.table); library(MMLB); library(MASS); library(parallel)
                  library(dplyr); library(fastGHQuad)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
for (m in 1:5) sys.source(file.path(SRC, sprintf("generation_model%d.R", m)), envir=globalenv())
sys.source(file.path(SRC,"generate_data.R"), envir=globalenv())
sys.source(file.path(SRC,"make_all_data.R"), envir=globalenv())

gen_gby <- function(m, N, T, i, label) {
  set.seed(20260805L + 1000L*m + i)
  as.data.frame(make_all_data(N, T, i, as.character(m), label))
}
fit_m3 <- function(d) suppressWarnings(mm(binY ~ time*trt, lv.formula=~time, t.formula=~1,
                                          data=d, id=id, q=20))

cat("===== VERIFICA 1: single fit GenBinaryY M3, N=60, T=6, replica 1 =====\n")
for (lb in c("mmlb","manual")) {
  d <- gen_gby(3, 60, 6, 1, lb)
  mo <- tryCatch(fit_m3(d), error=function(e) NULL)
  if (is.null(mo)) { cat("label", lb, ": fit fallito\n"); next }
  tb <- summary(mo)$mean.table
  k  <- which(rownames(tb)=="time:trt")
  cat(sprintf("\n-- label = %s (nobs=%d, eventi=%d)\n", lb, nrow(d), sum(d$binY)))
  cat("colnames(mean.table):", paste(colnames(tb), collapse=" | "), "\n")
  cat("rownames(mean.table):", paste(rownames(tb), collapse=" | "), "\n")
  cat(sprintf("time:trt  est            = %.10f\n", tb[k,1]))
  cat(sprintf("time:trt  mean.table SE  = %.10f\n", tb[k,2]))
  cat(sprintf("time:trt  sqrt(mod.cov)  = %.10f\n", sqrt(diag(mo$mod.cov))[k]))
  cat(sprintf("time:trt  sqrt(rob.cov)  = %.10f\n", sqrt(diag(mo$rob.cov))[k]))
  cat(sprintf("|tab - mod| = %.3e ; |tab - rob| = %.3e\n",
      abs(tb[k,2]-sqrt(diag(mo$mod.cov))[k]), abs(tb[k,2]-sqrt(diag(mo$rob.cov))[k])))
  cat("control[7] (cond.like/... settimo elemento) =", mo$control[7],
      " nome:", names(mo$control)[7], "\n")
}

cat("\n===== VERIFICA 2 (M-30): calibrazione a N=60, T=6, 400 repliche =====\n")
one <- function(i, lb) {
  d  <- tryCatch(gen_gby(3, 60, 6, i, lb), error=function(e) NULL)
  if (is.null(d)) return(NULL)
  mo <- tryCatch(fit_m3(d), error=function(e) NULL)
  if (is.null(mo)) return(NULL)
  tb <- tryCatch(summary(mo)$mean.table, error=function(e) NULL)
  if (is.null(tb) || !"time:trt" %in% rownames(tb)) return(NULL)
  k <- which(rownames(tb)=="time:trt")
  data.table(rep=i, label=lb, est=tb[k,1], se_tab=tb[k,2],
             se_mod=sqrt(diag(mo$mod.cov))[k], se_rob=sqrt(diag(mo$rob.cov))[k])
}
for (lb in c("mmlb","manual")) {
  rr <- rbindlist(Filter(Negate(is.null), mclapply(1:400, function(i) one(i, lb), mc.cores=40)))
  U <- rr[is.finite(est) & is.finite(se_mod) & abs(est)<10 & se_mod>0 & se_mod<10]
  cat(sprintf("\n-- label=%s : fit presenti %d/400 ; utilizzabili %d\n", lb, nrow(rr), nrow(U)))
  cat(sprintf("   SE media model-based = %.5f ; SD empirica = %.5f ; rapporto = %.4f\n",
      mean(U$se_mod), sd(U$est), mean(U$se_mod)/sd(U$est)))
  Ur <- U[is.finite(se_rob) & se_rob>0 & se_rob<10]
  cat(sprintf("   SE media robusta     = %.5f (n=%d) ; rapporto = %.4f\n",
      mean(Ur$se_rob), nrow(Ur), mean(Ur$se_rob)/sd(Ur$est)))
  cat(sprintf("   rapporto medio se_rob/se_mod = %.4f\n", mean(Ur$se_rob/Ur$se_mod)))
  fwrite(rr, path.expand(sprintf("~/MMM/rev2/VER_mmmrob_n60_%s.csv", lb)))
}
cat("\nDONE_VER_MMMROB\n")
