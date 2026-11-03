# REV2M_l.R -- remaining Methods-adjacent numbers: D-07's 44,169; §3.4's 8.51/9.17;
# the LRT/Wald cell-count arithmetic; §2.4's exclusion-rate wording.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt"]; if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  x[, .(dgp,method,model,N,ntime,label,sim_id,converged,estimate,se,p)]
}
L <- mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)
D <- rbindlist(L[!sapply(L,is.null)], fill=TRUE)
cat("total replicate rows loaded:", nrow(D), "\n")

cat("\n########## D-07: LRT rows a Wald SE filter would discard\n")
lr <- D[method=="glmer_lrt" & converged & is.finite(p)]
cat(sprintf("LRT usable rows under the per-method filter: %d\n", nrow(lr)))
drop <- lr[!(is.finite(se) & se>0 & se<10 & is.finite(estimate) & abs(estimate)<10)]
cat(sprintf("of those, rows a Wald SE/estimate filter would drop: %d\n", nrow(drop)))
print(as.data.frame(drop[, .(dropped=.N), by=dgp]), row.names=FALSE)
cat(sprintf("level with the per-method filter: %.2f%% | with a Wald filter applied to LRT rows: %.2f%%\n",
  100*mean(lr[label=="manual"]$p<0.05),
  100*mean(lr[label=="manual"][is.finite(se) & se>0 & se<10 & abs(estimate)<10]$p<0.05)))
cat(sprintf("difference: %.2f pp (claim: 10.92 -> 9.32, i.e. 1.60 pp; manuscript says 1.6 pp)\n",
  100*mean(lr[label=="manual"]$p<0.05) -
  100*mean(lr[label=="manual"][is.finite(se) & se>0 & se<10 & abs(estimate)<10]$p<0.05)))
cat("\nsame per mechanism:\n")
z <- lr[label=="manual", .(all=100*mean(p<0.05),
   waldfiltered=100*mean(p[is.finite(se) & se>0 & se<10 & abs(estimate)<10]<0.05),
   n_all=.N, n_wald=sum(is.finite(se) & se>0 & se<10 & abs(estimate)<10)), by=dgp]
z[, drop_pp := round(all-waldfiltered,2)]
print(as.data.frame(z[, lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)]), row.names=FALSE)

cat("\n########## §3.4: shrinkage vs GLMM Wald rejection rate on GenBinaryY\n")
u <- D[converged & is.finite(estimate) & is.finite(se) & is.finite(p) & abs(estimate)<10 & se>0 & se<10]
w <- u[label=="manual" & method %in% c("glmer","shrink")]
cat("own-set worst-case / pooled level:\n")
print(as.data.frame(w[, .(n=.N, level=round(100*mean(p<0.05),2)), by=.(dgp,method)][order(dgp,method)]),
      row.names=FALSE)
cat("\npaired on replicates where BOTH converged:\n")
pw <- dcast(w[, .(dgp,model,N,ntime,sim_id,method,p)], dgp+model+N+ntime+sim_id~method, value.var="p")
pw <- pw[!is.na(glmer) & !is.na(shrink)]
print(as.data.frame(pw[, .(n=.N, glmer=round(100*mean(glmer<0.05),2), shrink=round(100*mean(shrink<0.05),2),
   discordant=sum((glmer<0.05)!=(shrink<0.05))), by=dgp]), row.names=FALSE)
cat("\n########## exclusion of |estimate|>=10 : rate per method and mechanism\n")
cv <- D[converged==TRUE]
ex <- cv[, .(converged=.N, big=sum(is.finite(estimate) & abs(estimate)>=10),
             no_se=sum(!is.finite(se)), se_big=sum(is.finite(se) & se>=10)), by=.(method,dgp)]
ex[, `:=`(pct_big=round(100*big/converged,2), pct_nose=round(100*no_se/converged,2),
          pct_sebig=round(100*se_big/converged,3))]
print(as.data.frame(ex[order(-pct_big)][1:15]), row.names=FALSE)
cat("\n########## LRT/Wald specification mismatch: cell arithmetic\n")
cat(sprintf("scenarios with a random slope in the Wald arm: M3,M4,M5 = 3 of 5\n"))
cat(sprintf("design points affected: 3 scenarios x 3 mechanisms x 16 N x 6 T = %d of %d (%.1f%%)\n",
  3*3*16*6, 5*3*16*6, 100*(3*3*16*6)/(5*3*16*6)))
cat(sprintf("cells (design point x hypothesis) affected: %d of %d\n", 2*3*3*16*6, 2*5*3*16*6))
cat(sprintf("replicate rows of glmer_lrt in the affected scenarios: %d of %d (%.1f%%)\n",
  nrow(lr[model %in% 3:5]), nrow(lr), 100*nrow(lr[model %in% 3:5])/nrow(lr)))
cat("\n########## LRT usable share at T=2 by scenario group (the 98.8% figure)\n")
att <- D[, .(n_attempt=uniqueN(sim_id)), by=.(dgp,model,N,ntime,label)]
lrc <- lr[, .(usable=.N), by=.(dgp,model,N,ntime,label)]
mg <- merge(att, lrc, by=c("dgp","model","N","ntime","label"), all.x=TRUE)
mg[is.na(usable), usable := 0L]
mg[, re := ifelse(model %in% c(1,2), "RI", "RS")]
print(as.data.frame(mg[ntime==2, .(attempts=sum(n_attempt), usable=sum(usable),
   pct=round(100*sum(usable)/sum(n_attempt),2)), by=re]), row.names=FALSE)
print(as.data.frame(mg[, .(pct=round(100*sum(usable)/sum(n_attempt),2)), by=.(re,ntime)][order(re,ntime)]),
      row.names=FALSE)
cat("\nDONE_L\n")
