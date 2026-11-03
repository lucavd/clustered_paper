.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(MMLB); library(lme4); library(geepack)
                  library(geesmv); library(dplyr); library(tibble); library(MASS); library(fastGHQuad)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")

# (A) esattamente come il worker di run_full.R: sys.source in globalenv
sp <- function(f, dir=RD) sys.source(file.path(dir, f), envir=globalenv())
for (f in c("analysis_model1.R","analysis_model3.R","do_analysis.R")) sp(f)
for (f in c("generation_model1.R","generation_model3.R","generate_data.R","make_all_data.R")) sp(f)

ch <- as.data.table(readRDS(path.expand("~/MMM/results_full/GenBinaryY_m1_N10_T4_manual_c01.rds")))
cm <- ch[term=="time:trt" & method=="mmm"][order(sim_id)]
cat("chunk: righe mmm", nrow(cm), " sim_id presenti:", paste(head(cm$sim_id,8),collapse=","), "\n")

cat("\n-- (A) globalenv, seed 20260805+1000*model+i, model INTERO (come run_full) --\n")
for (i in 1:5) {
  set.seed(20260805L + 1000L*1L + i)
  d <- tryCatch(make_all_data(10, 4, i, 1L, "manual"), error=function(e) NULL)
  if (is.null(d)) { cat(sprintf("i=%d generazione fallita\n", i)); next }
  d <- as.data.frame(d)
  mo <- tryCatch(suppressWarnings(mm(binY ~ time*trt, lv.formula=~1, t.formula=~1, data=d, id=id, q=20)),
                 error=function(e) NULL)
  e <- if (is.null(mo)) NA_real_ else summary(mo)$mean.table["time:trt",1]
  cz <- cm[sim_id==i]$estimate
  cat(sprintf("i=%d  nobs=%3d eventi=%2d  est_A=%s  est_chunk=%s\n", i, nrow(d), sum(d$binY),
              formatC(e, format="f", digits=8),
              if (!length(cz)) "assente" else formatC(cz, format="f", digits=8)))
}
cat("\n-- (B) come sopra ma model come CARATTERE '1' --\n")
for (i in 1:3) {
  set.seed(20260805L + 1000L*1L + i)
  d <- as.data.frame(make_all_data(10, 4, i, "1", "manual"))
  cat(sprintf("i=%d nobs=%d eventi=%d  primi binY: %s\n", i, nrow(d), sum(d$binY),
              paste(head(d$binY,10), collapse="")))
}
cat("\n-- (C) la stessa generazione ripetuta due volte: determinismo --\n")
set.seed(20260805L + 1001L); d1 <- as.data.frame(make_all_data(10, 4, 1, 1L, "manual"))
set.seed(20260805L + 1001L); d2 <- as.data.frame(make_all_data(10, 4, 1, 1L, "manual"))
cat("identiche:", identical(d1$binY, d2$binY), " eventi:", sum(d1$binY), sum(d2$binY), "\n")

cat("\n-- (D) la generazione consuma RNG dopo la seed? conteggio estrazioni --\n")
set.seed(1); a <- runif(1)
set.seed(20260805L + 1001L); invisible(make_all_data(10, 4, 1, 1L, "manual"))
cat("stato RNG dopo generazione, prossimo runif:", runif(1), "\n")

cat("\n-- (E) semi alternativi contro il chunk (i=1..3) --\n")
cands <- list("20260805+1000*m+i" = function(m,i) 20260805L + 1000L*m + i,
              "20260805+i"        = function(m,i) 20260805L + i,
              "i"                 = function(m,i) i,
              "20260805*m+i"      = function(m,i) 20260805L*m + i,
              "1000*m+i"          = function(m,i) 1000L*m + i)
for (nm in names(cands)) {
  out <- character(0)
  for (i in 1:3) {
    set.seed(cands[[nm]](1L, i))
    d <- as.data.frame(make_all_data(10, 4, i, 1L, "manual"))
    mo <- tryCatch(suppressWarnings(mm(binY ~ time*trt, lv.formula=~1, t.formula=~1, data=d, id=id, q=20)),
                   error=function(e) NULL)
    e <- if (is.null(mo)) NA_real_ else summary(mo)$mean.table["time:trt",1]
    out <- c(out, sprintf("%.6f", e))
  }
  cat(sprintf("%-20s : %s\n", nm, paste(out, collapse="  ")))
}
cat("chunk (i=1..3)        : ", paste(sprintf("%.6f", cm[sim_id %in% 1:3]$estimate), collapse="  "), "\n")
cat("\n-- (F) il chunk contiene sim_id 1..100? --\n")
cat("range sim_id chunk:", range(ch$sim_id), " unici:", uniqueN(ch$sim_id), "\n")
cat("DONE_SEEDCHK\n")
