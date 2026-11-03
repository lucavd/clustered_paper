
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
fs <- list.files(RES, pattern="^GenBinaryY_.*manual_c[0-9]+[.]rds$", full.names=TRUE)
k <- sub("_c[0-9]+[.]rds$","",basename(fs)); fs <- fs[k %in% names(which(table(k)==10))]
cat("H0 GenBinaryY files:", length(fs), "\n")
D <- rbindlist(mclapply(fs, function(x){ z <- as.data.table(readRDS(x))
  z[term=="time:trt", .(estimate,se,stat,p,method,N,ntime,model,converged)] }, mc.cores=40))
D[, model := as.integer(model)]
H0 <- rbind(D[method %in% WALD & converged & is.finite(estimate) & is.finite(se) &
              abs(estimate)<10 & se>0 & se<10 & is.finite(p)],
            D[method=="glmer_lrt" & converged & is.finite(p)])

cat("\n### is p a monotone function of |estimate/se| within method? (H0 rows)\n")
H0[method!="glmer_lrt", z := abs(estimate/se)]
print(H0[method!="glmer_lrt", .(n=.N, spearman_p_vs_absz=round(cor(p, z, method="spearman"),6),
        max_abs_dev = round(max(abs(p - 2*pnorm(-z))), 8)), by=method][order(method)])

cat("\n### ties / degeneracy in the p vector, per method (H0 GenBinaryY)\n")
print(H0[, .(n=.N, n_distinct=uniqueN(p), n_p_eq_0=sum(p==0), n_p_eq_1=sum(p==1),
             min_p=format(min(p),digits=3)), by=method][order(method)])

cat("\n### per-cell empirical size when q05 is applied back to the same H0 rows\n")
CS <- H0[, { q <- quantile(p, 0.05, na.rm=TRUE, type=7)
             .(nH0=.N, q05=q, size=100*mean(p<q), size_le=100*mean(p<=q),
               ties_at_q=sum(p==q)) }, by=.(method,model,N,ntime)]
print(CS[, .(ncell=.N, size_pooled_rowwt=round(sum(size*nH0/100)/sum(nH0)*100,4),
             size_cellmean=round(mean(size),4), min_cell=round(min(size),3),
             max_cell=round(max(size),3), cells_below_4.9=sum(size<4.9),
             cells_with_ties_at_q=sum(ties_at_q>1)), by=method][order(method)])
cat("\n### glmer/shrink: cells with size < 4.9, by number of usable H0 rows\n")
print(CS[method %in% c("glmer","shrink") & size<4.9,
         .(ncell=.N, median_nH0=median(nH0), min_nH0=min(nH0), max_nH0=max(nH0),
           median_size=round(median(size),3), median_ties=median(ties_at_q)), by=method])
cat("\n### size vs nH0 decile, glmer\n")
G <- CS[method=="glmer"][order(nH0)]
G[, dec := cut(nH0, breaks=quantile(nH0, 0:5/5), include.lowest=TRUE)]
print(G[, .(ncell=.N, nH0_med=median(nH0), size_mean=round(mean(size),3),
            ties_med=median(ties_at_q)), by=dec])
cat("\n### using p <= q05 instead of p < q05 (tie handling), pooled\n")
print(CS[, .(size_lt=round(sum(size*nH0)/sum(nH0),4),
             size_le=round(sum(size_le*nH0)/sum(nH0),4)), by=method][order(method)])
fwrite(CS, file.path(OUT,"R2_size_percell_GenBinaryY.csv"))
cat("\nDONE_D\n")
