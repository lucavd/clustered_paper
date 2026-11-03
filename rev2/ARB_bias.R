
# Arbitrato indipendente del bias condizionale del GLMM (sezione 3.13), dai CHUNK,
# non dai CSV intermedi. Quattro rilievi danno 0,0460 / 0,0453 / 0,0525 / 0,0798.
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")

# bersagli in circolazione
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)   # marginale per scenario
COND_FIX <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232) # condizionale fisso
AR1 <- fread(path.expand("~/MMM/rev2/R2_ar1_cond_true.csv"))     # cond. calibrato AR1 per (model,T)
setnames(AR1, "t_cond_true", "tc")

fs <- list.files(RES, pattern="_mmlb_c[0-9]+[.]rds$")
cat("chunk H1 totali:", length(fs), "\n")
one <- function(f) {
  x <- try(as.data.table(readRDS(file.path(RES,f))), silent=TRUE)
  if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method=="glmer"]
  if (!nrow(x)) return(NULL)
  x[, converged := as.logical(converged)]; x[is.na(converged), converged := FALSE]
  u <- x[converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
         abs(estimate) < 10 & se > 0 & se < 10]
  if (!nrow(u)) return(NULL)
  u[, .(dgp, model=as.integer(model), N, ntime, estimate)]
}
L <- mclapply(fs, one, mc.cores=40)
D <- rbindlist(Filter(is.data.table, L))
cat("repliche glmer utilizzabili:", nrow(D), " celle:", uniqueN(D[,.(dgp,model,N,ntime)]), "\n")

D[, t_marg := MARG[as.character(model)]]
D[, t_cfix := COND_FIX[as.character(model)]]
D <- merge(D, AR1[, .(model, ntime, tc)], by=c("model","ntime"), all.x=TRUE)
# convenzione "calibrata": AR1 usa il valore per (model,T); gli altri due il fisso
D[, t_ccal := ifelse(dgp=="GLMMAR1" | dgp=="GLMM-AR1", tc, t_cfix)]
cat("AR1 senza calibrazione:", D[(dgp %like% "AR1") & is.na(tc), .N], "\n")

f <- function(v, lab) {
  # pesato per repliche, e media per cella, del bias SEGNATO e ASSOLUTO
  w  <- D[, .(bias = mean(estimate - get(v)))]
  pc <- D[, .(b = mean(estimate - get(v))), by=.(dgp,model,N,ntime)][, .(bias = mean(b), abs = mean(abs(b)))]
  wa <- D[, .(b = mean(estimate - get(v))), by=.(dgp,model,N,ntime)]
  cat(sprintf("%-26s segnato pesato %.6f | segnato per cella %.6f | |bias| per cella %.6f  (celle %d)\n",
              lab, w$bias, pc$bias, pc$abs, nrow(wa)))
}
cat("\n=== BIAS DEL GLMM, n =", nrow(D), "repliche ===\n")
f("t_marg","marginale per scenario")
f("t_cfix","condizionale FISSO 0.232/0.058")
f("t_ccal","condizionale CALIBRATO per AR1")

# riduzione percentuale del |bias| passando dal marginale al condizionale, per cella
red <- function(v) {
  a <- D[, .(bm = mean(estimate - t_marg), bc = mean(estimate - get(v))), by=.(dgp,model,N,ntime)]
  100*(mean(abs(a$bm)) - mean(abs(a$bc)))/mean(abs(a$bm))
}
cat(sprintf("\nriduzione del |bias| per cella, condizionale FISSO:     %.4f%%\n", red("t_cfix")))
cat(sprintf("riduzione del |bias| per cella, condizionale CALIBRATO: %.4f%%\n", red("t_ccal")))
cat("\nper meccanismo (condizionale calibrato, segnato pesato):\n")
print(D[, .(n=.N, bias=round(mean(estimate - t_ccal),6)), by=dgp])
