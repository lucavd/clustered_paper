
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)
CFIX <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232)
F1 <- fread(path.expand("~/MMM/tables/F_truth_ar1.csv"))[, .(tc = mean(cond_used)), by=.(model,ntime)]
fs <- list.files(RES, pattern="_mmlb_c[0-9]+[.]rds$")
one <- function(f) {
  x <- try(as.data.table(readRDS(file.path(RES,f))), silent=TRUE)
  if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method=="glmer"]; if (!nrow(x)) return(NULL)
  x[, converged := as.logical(converged)]; x[is.na(converged), converged := FALSE]
  u <- x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10]
  if (!nrow(u)) return(NULL)
  u[, .(dgp, model=as.integer(model), N, ntime, estimate)]
}
D <- rbindlist(Filter(is.data.table, mclapply(fs, one, mc.cores=40)))
D <- merge(D, F1, by=c("model","ntime"), all.x=TRUE)
D[, t_marg := MARG[as.character(model)]]
D[, t_cond := ifelse(dgp %like% "AR1", tc, CFIX[as.character(model)])]
pc <- D[, .(bm=mean(estimate-t_marg), bc=mean(estimate-t_cond)), by=.(dgp,model,N,ntime)]
cat(sprintf("FONTE F_truth_ar1.csv (tracciata nel repo), n=%d repliche, %d celle\n", nrow(D), nrow(pc)))
cat(sprintf("  bias condizionale: pesato %.6f | per cella %.6f | |bias| per cella %.6f\n",
            D[, mean(estimate-t_cond)], mean(pc$bc), mean(abs(pc$bc))))
cat(sprintf("  bias marginale   : pesato %.6f | per cella %.6f | |bias| per cella %.6f\n",
            D[, mean(estimate-t_marg)], mean(pc$bm), mean(abs(pc$bm))))
cat(sprintf("  riduzione del |bias| per cella: %.4f%%\n",
            100*(mean(abs(pc$bm))-mean(abs(pc$bc)))/mean(abs(pc$bm))))
