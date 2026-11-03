# REV2M_a6.R -- §3.8 / §3.9 tables under BOTH denominators, and the LRT level range.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table)})
R <- fread(path.expand("~/MMM/rev2/REV_cells.csv")); R[, model := as.integer(model)]
D <- unique(R[set=="own", .(dgp,model,N,ntime,label,method,n_attempt,n_usable)])
METH <- c("glmer","shrink","glmer_lrt","mmm","gee_fg")
# full cell frame: every (dgp,model,N,ntime,label) x method, missing = 0 usable, 1000 attempts
cells <- unique(R[, .(dgp,model,N,ntime,label)])
full <- CJ(k=1:nrow(cells), method=METH)[, cbind(cells[k], method=method)]
D2 <- merge(full, D[, .(dgp,model,N,ntime,label,method,n_usable,n_attempt)],
            by=c("dgp","model","N","ntime","label","method"), all.x=TRUE)
D2[is.na(n_usable), n_usable := 0L]; D2[is.na(n_attempt), n_attempt := 1000L]
D2[, re := ifelse(model %in% c(1,2), "RI", "RS")]
cat("cells in the full frame per method:", nrow(D2)/length(METH), "\n\n")
cat("===== §3.8 usable %, denominator = 1000 x cells (a never-fitting cell counts 0/1000) =====\n")
t1 <- D2[, .(pct=round(100*sum(n_usable)/(1000*.N),2)), by=.(re,ntime,method)]
print(as.data.frame(dcast(t1, re+ntime~method, value.var="pct")), row.names=FALSE)
cat("\n===== same, denominator = replicates actually present (n_attempt) =====\n")
t2 <- D2[, .(pct=round(100*sum(n_usable)/sum(n_attempt),2)), by=.(re,ntime,method)]
print(as.data.frame(dcast(t2, re+ntime~method, value.var="pct")), row.names=FALSE)
cat("\n===== pooled over all 5 scenarios at T=2 =====\n")
print(as.data.frame(D2[ntime==2, .(usable=sum(n_usable), denom1000=1000*.N,
  pct_1000=round(100*sum(n_usable)/(1000*.N),2),
  pct_present=round(100*sum(n_usable)/sum(n_attempt),2)), by=method]), row.names=FALSE)
cat("\n===== RS, T=2, per arm, raw counts against 144,000 =====\n")
print(as.data.frame(D2[re=="RS" & ntime==2, .(cells=.N, denom=1000*.N,
   usable=sum(n_usable)), by=.(method,label)][order(method,label)]), row.names=FALSE)
cat("\n===== §3.9: N<=20, 270 design points per method, cells never fitted counted at 0 =====\n")
ALL8 <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink","glmer_lrt")
full8 <- CJ(k=1:nrow(cells), method=ALL8)[, cbind(cells[k], method=method)]
D8 <- merge(full8, D, by=c("dgp","model","N","ntime","label","method"), all.x=TRUE)
D8[is.na(n_usable), n_usable := 0L]; D8[is.na(n_attempt), n_attempt := 1000L]
K <- unique(R[set=="own", .(dgp,model,N,ntime,label,method,n,k_rej)])
D8 <- merge(D8, K, by=c("dgp","model","N","ntime","label","method"), all.x=TRUE)
D8[is.na(n), n := 0L]; D8[is.na(k_rej), k_rej := 0L]
s <- D8[N<=20 & label=="manual"]
out <- s[, .(design_points=uniqueN(paste(dgp,model,N,ntime)),
             usable_pct = round(100*mean(n_usable/1000),2),
             val_pct = round(100*mean((n_usable/1000)*ifelse(n>0, 1-k_rej/n, 0)),2),
             dp_never = uniqueN(paste(dgp,model,N,ntime)[n_usable==0])), by=method][order(-usable_pct)]
print(as.data.frame(out), row.names=FALSE)
cat("\n===== the same restricted to cells each method attempts =====\n")
out2 <- s[n_usable>0, .(design_points=uniqueN(paste(dgp,model,N,ntime)),
   usable_pct=round(100*mean(n_usable/1000),2)), by=method][order(-usable_pct)]
print(as.data.frame(out2), row.names=FALSE)
cat("\n===== LRT level at N<=20: replicate-weighted, and mean over design points =====\n")
l1 <- s[method=="glmer_lrt", .(rep_weighted=round(100*sum(k_rej)/sum(n),2),
        dp_mean=round(100*mean(ifelse(n>0,k_rej/n,NA), na.rm=TRUE),2)), by=dgp]
print(as.data.frame(l1), row.names=FALSE)
cat(sprintf("range replicate-weighted: %.2f - %.2f ; range design-point mean: %.2f - %.2f\n",
  min(l1$rep_weighted), max(l1$rep_weighted), min(l1$dp_mean), max(l1$dp_mean)))
cat("\n===== corrected sandwich usable %, N<=20 and whole grid =====\n")
cs <- D8[label=="manual" & method %in% c("gee_fg","gee_kc","gee_md","glmer_lrt")]
print(as.data.frame(cs[, .(N_le20=round(100*mean(n_usable[N<=20]/1000),2),
                           whole_grid=round(100*mean(n_usable/1000),2)), by=method]), row.names=FALSE)
cat("\nDONE_A6\n")
